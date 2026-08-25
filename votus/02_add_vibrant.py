#!/usr/bin/env python3
"""
votus/02_add_vibrant.py — VIBRANT lifestyle and genome quality per vOTU.

Rebuild of 02_add_vibrant.py, with two things fixed:

  1. FRAGMENT NAMES. VIBRANT reports excised prophages as
     <scaffold>_fragment_N, which never matched the catalogue's vOTU IDs.
     Exactly the problem that lost 1,290 ORFs in the gene branch. Here the
     suffix is stripped so every prediction lands on its parent contig.

  2. DUPLICATES RESOLVED ON PURPOSE. The original did
         drop_duplicates(subset=['vOTU_ID'], keep='first')
     with the comment "safety patch". But a contig carrying two prophages
     legitimately produces two VIBRANT rows, and 'first' picks one at
     random — possibly assigning the wrong lifestyle. Now:
         quality  -> the BEST of the rows
         lifestyle-> lysogenic wins over lytic, because a prophage
                     prediction is positive evidence and its absence is not
     and both the number of collapsed rows and the number of conflicting
     lifestyles are counted and printed.

  3. Also adds Vibrant_predicted (Yes/No): whether VIBRANT called this
     contig viral at all. That is one of the five independent lines of
     evidence used by the catalogue filter in step 06, and it was being
     inferred indirectly from the Quality string before.

Output: votu_02_vibrant.tsv
"""

import os
import re
import sys

import pandas as pd

from config import F_VOTU_01, F_VOTU_02, F_VIB_QUALITY, require, Counter_

FRAG_RE = re.compile(r"_fragment_\d+$")

# Best first. Anything unrecognised sorts last and is counted.
QUALITY_RANK = [
    "complete circular",
    "high quality draft",
    "medium quality draft",
    "low quality draft",
]

cnt = Counter_("COUNTERS — nothing is dropped in silence")


def q_rank(q):
    s = str(q).strip().lower()
    return QUALITY_RANK.index(s) if s in QUALITY_RANK else len(QUALITY_RANK)


def main():
    require(F_VOTU_01, F_VIB_QUALITY)
    print("=" * 70)
    print(" vOTU LEVEL 2 — VIBRANT lifestyle and quality")
    print("=" * 70)

    df = pd.read_csv(F_VOTU_01, sep="\t", dtype=str)
    n0 = len(df)
    print(f"[1] base table: {n0} vOTUs")

    # ---------------------------------------------------------------- 2
    print("[2] reading VIBRANT genome quality ...")
    vq = pd.read_csv(F_VIB_QUALITY, sep="\t", dtype=str)
    need = ["scaffold", "type", "Quality"]
    missing = [c for c in need if c not in vq.columns]
    if missing:
        raise SystemExit(f"[SCHEMA ERROR] VIBRANT quality lacks {missing}\n"
                         f"  found: {list(vq.columns)}")
    vq = vq[need].copy()
    print(f"    VIBRANT predictions: {len(vq)}")

    print("    observed values:")
    for c in ("type", "Quality"):
        vals = vq[c].fillna("(empty)").value_counts()
        print(f"      {c}: " + ", ".join(f"{k} ({v})" for k, v in vals.items()))
    unknown_q = sorted({q for q in vq["Quality"].fillna("")
                        if str(q).strip().lower() not in QUALITY_RANK})
    if unknown_q:
        cnt.add(f"Quality values not in QUALITY_RANK: {unknown_q}")

    # ---------------------------------------------------------------- 3
    print("[3] mapping excised-prophage names back to the parent contig ...")
    vq["is_fragment"] = vq["scaffold"].str.contains("_fragment_", na=False)
    n_frag = int(vq["is_fragment"].sum())
    vq["vOTU_ID"] = vq["scaffold"].str.replace(FRAG_RE, "", regex=True)
    print(f"    predictions on excised prophages: {n_frag}")
    cnt.add("VIBRANT predictions named _fragment_N (remapped to parent)", n_frag)

    # ---------------------------------------------------------------- 4
    print("[4] resolving contigs with more than one prediction ...")
    vq["_qrank"] = vq["Quality"].apply(q_rank)
    vq["_lyso"] = vq["type"].str.strip().str.lower().eq("lysogenic")

    grp = vq.groupby("vOTU_ID", sort=True)
    n_multi = int((grp.size() > 1).sum())
    print(f"    contigs with >1 VIBRANT prediction: {n_multi}")
    cnt.add("contigs with >1 VIBRANT prediction (collapsed deliberately)",
            n_multi)

    conflict = int(grp["_lyso"].nunique().gt(1).sum())
    if conflict:
        print(f"    of those, lifestyle CONFLICT (lytic + lysogenic): {conflict}")
        print("      -> resolved as lysogenic: a prophage prediction is "
              "positive evidence, its absence is not")
        cnt.add("lifestyle conflicts resolved as lysogenic", conflict)

    res = pd.DataFrame({
        "Vibrant_Quality": grp.apply(
            lambda g: g.sort_values("_qrank")["Quality"].iloc[0]),
        "Vibrant_Lifestyle": grp["_lyso"].any().map(
            {True: "lysogenic", False: "lytic"}),
        "Vibrant_n_predictions": grp.size(),
        "Vibrant_from_fragment": grp["is_fragment"].any().map(
            {True: "Yes", False: "No"}),
    }).reset_index()

    # ---------------------------------------------------------------- 5
    print("[5] merging ...")
    df = df.merge(res, on="vOTU_ID", how="left")
    if len(df) != n0:
        print("[FATAL] the merge changed the row count — duplicate vOTU_IDs")
        sys.exit(1)

    df["Vibrant_predicted"] = df["Vibrant_Quality"].notna().map(
        {True: "Yes", False: "No"})
    n_pred = int((df["Vibrant_predicted"] == "Yes").sum())
    df["Vibrant_Quality"] = df["Vibrant_Quality"].fillna("Not-predicted")
    df["Vibrant_Lifestyle"] = df["Vibrant_Lifestyle"].fillna("Not-predicted")
    df["Vibrant_n_predictions"] = df["Vibrant_n_predictions"].fillna(0).astype(int)
    df["Vibrant_from_fragment"] = df["Vibrant_from_fragment"].fillna("No")

    orphan = set(res["vOTU_ID"]) - set(df["vOTU_ID"])
    if orphan:
        cnt.add("VIBRANT predictions whose contig is NOT in the catalogue",
                len(orphan))
        print(f"    [warn] {len(orphan)} predictions on contigs absent from "
              f"the catalogue, e.g. {sorted(orphan)[:2]}")

    df.to_csv(F_VOTU_02, sep="\t", index=False)

    # ---------------------------------------------------------------- report
    print("\n" + "=" * 70)
    print(" VERIFICATION")
    print("=" * 70)
    print(f" vOTUs in / out                 : {n0} / {len(df)}")
    print(f" predicted viral by VIBRANT     : {n_pred}  "
          f"({100.0*n_pred/n0:.1f}%)")
    print(f" NOT predicted by VIBRANT       : {n0-n_pred}")
    print("   (expected: the catalogue also admits contigs via CheckV,")
    print("    PhaBOX, CT3 and provirus evidence — VIBRANT is one of five)")

    print("\n VIBRANT lifestyle:")
    for k, v in df["Vibrant_Lifestyle"].value_counts().items():
        print(f"   {v:>8}  {k}")
    print("\n VIBRANT quality:")
    for k, v in df["Vibrant_Quality"].value_counts().items():
        print(f"   {v:>8}  {k}")

    print("\n lysogenic by culture group (VIBRANT only, no filtering yet):")
    sub = df[df["Vibrant_Lifestyle"] != "Not-predicted"]
    if len(sub):
        tab = pd.crosstab(sub["Fraction_group"], sub["Vibrant_Lifestyle"])
        for grp_name, row in tab.iterrows():
            tot = row.sum()
            lyso = row.get("lysogenic", 0)
            print(f"   {grp_name:<22} {lyso:>6} / {tot:<6} "
                  f"({100.0*lyso/tot:.1f}% lysogenic)")
        print("   [!] descriptive only — LIS is mitomycin-induced, so this")
        print("       contrast is confounded until CB22 vs GLI is checked")

    cnt.report(total=n0)
    print("\n" + "=" * 70)
    print(f" Output: {F_VOTU_02}")


if __name__ == "__main__":
    main()
