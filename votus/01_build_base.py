#!/usr/bin/env python3
"""
genes/01_build_gene_table.py — Build LEVEL 2 of the pipeline: one row per ORF.

The spine is the Prodigal gene call over the full vOTU catalogue 
(establishing ONE coordinate frame). VIBRANT annotations are layered 
on top as distinct, independent columns to ensure annotation fields 
(e.g., KO, Pfam, VOG) are never conflated during downstream parsing.

Prophage Fragment Re-mapping:
VIBRANT renames ORFs on excised prophages (e.g., <scaffold>_fragment_F_N).
Since the N index already represents the true ORF number relative to the 
parent contig (verified against the prophage coordinate file where 
min(N) == 'protein start' and max(N) == 'protein stop'), this script 
algorithmically remaps fragment ORFs back to their parent ID (<scaffold>_N). 
This guarantees seamless identity tracking across the pipeline.

Output: 01_gene_table_L2.tsv
"""

import os
import re
import sys

import pandas as pd

from config import (
    F_CATALOGUE_FNA, F_CLUSTER_TSV, F_CHECKV, F_VOTU_01,
    FRACTIONS, ZONES, require, Counter_,
)

NAME_RE = re.compile(r"MCHA-([A-Z]{2})-\d+-([A-Z0-9]+)")
PRO_RE = re.compile(r"_pro_?\d*$")
MAG_RE = re.compile(r"(MCHA-[A-Z]{2}-\d+-[A-Z0-9]+_bin\.\d+)")

cnt = Counter_("COUNTERS — nothing is dropped in silence")


def main():
    require(F_CATALOGUE_FNA, F_CLUSTER_TSV, F_CHECKV)
    print("=" * 70)
    print(" vOTU LEVEL 1 — identity, sampling design, cluster, CheckV")
    print("=" * 70)

    # ---------------------------------------------------------------- 1
    print("[1] reading vOTU IDs from the catalogue FASTA ...")
    ids = []
    with open(F_CATALOGUE_FNA) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0].strip())
    df = pd.DataFrame({"vOTU_ID": ids})
    n0 = len(df)
    print(f"    vOTUs: {n0}")
    if df["vOTU_ID"].duplicated().any():
        dup = int(df["vOTU_ID"].duplicated().sum())
        cnt.add("duplicate IDs in the FASTA (kept, but CHECK THIS)", dup)

    # ---------------------------------------------------------------- 2
    print("[2] sampling design from the sequence names ...")

    def design(vid):
        m = NAME_RE.search(vid)
        if not m:
            cnt.add("name did not match MCHA-ZZ-n-FRACTION")
            return pd.Series(["Unknown", "Unknown", "Unknown",
                              "Unknown", "Unknown", "Unknown"])
        zone, frac = m.group(1), m.group(2)
        if zone not in ZONES:
            cnt.add(f"zone not in config.ZONES: {zone}")
        f = FRACTIONS.get(frac)
        if f is None:
            cnt.add(f"fraction not in config.FRACTIONS: {frac}")
            return pd.Series([zone, frac, "Unknown",
                              "Unknown", "Unknown", "Unknown"])
        return pd.Series([zone, frac, f["group"], f["mda"],
                          f["mitomycin"], f["quantified_own"]])

    df[["Zone", "Fraction", "Fraction_group",
        "MDA", "Mitomycin", "Quantified_own"]] = df["vOTU_ID"].apply(design)

    print("    fractions found:")
    for frac, n in df["Fraction"].value_counts().items():
        f = FRACTIONS.get(frac, {})
        print(f"      {frac:<8} {n:>7}   {f.get('group','?'):<22}"
              f" MDA={f.get('mda','?'):<4} mitomycin={f.get('mitomycin','?')}")

    # ---------------------------------------------------------------- 3
    print("[3] tracing clusters for provirus and MAG origin ...")
    clu = pd.read_csv(F_CLUSTER_TSV, sep="\t", header=None,
                      names=["Rep", "Member"], dtype=str)
    clu["Member"] = clu["Member"].fillna("")

    rows = []
    for rep, g in clu.groupby("Rep", sort=True):
        members = g["Member"].tolist()
        is_pro = any(PRO_RE.search(m) for m in members)
        mags = sorted({m.group(1) for x in members
                       for m in [MAG_RE.search(x)] if m})   # sorted == stable
        rows.append({
            "vOTU_ID": rep,
            "Cluster_size": len(members),
            "Cluster_has_Provirus": "Yes" if is_pro else "No",
            "MAG_n": len(mags),
            "MAG_ids": ";".join(mags) if mags else "None",
        })
    clusters = pd.DataFrame(rows)
    print(f"    clusters: {len(clusters)}")

    df = df.merge(clusters, on="vOTU_ID", how="left")
    miss = int(df["Cluster_size"].isna().sum())
    if miss:
        cnt.add("vOTUs absent from the cluster file (singletons assumed)", miss)
    df["Cluster_size"] = df["Cluster_size"].fillna(1).astype(int)
    df["Cluster_has_Provirus"] = df["Cluster_has_Provirus"].fillna("No")
    df["MAG_n"] = df["MAG_n"].fillna(0).astype(int)
    df["MAG_ids"] = df["MAG_ids"].fillna("None")

    n_mag = int((df["MAG_n"] > 0).sum())
    n_multi = int((df["MAG_n"] > 1).sum())
    print(f"    vOTUs linked to >=1 MAG : {n_mag}")
    print(f"    vOTUs linked to >1 MAG  : {n_multi}  "
          f"(all listed; taxonomy resolved in step 04)")
    cnt.add("vOTUs linked to more than one MAG", n_multi)

    # ---------------------------------------------------------------- 4
    print("[4] adding CheckV ...")
    cv = pd.read_csv(F_CHECKV, sep="\t", dtype=str)
    need = ["contig_id", "contig_length", "checkv_quality",
            "completeness", "provirus"]
    missing = [c for c in need if c not in cv.columns]
    if missing:
        raise SystemExit(f"[SCHEMA ERROR] CheckV lacks columns: {missing}")
    cv = cv[need].rename(columns={"contig_id": "vOTU_ID",
                                  "provirus": "CheckV_Provirus"})
    cv = cv.drop_duplicates(subset=["vOTU_ID"], keep="first")

    df = df.merge(cv, on="vOTU_ID", how="left")
    miss_cv = int(df["checkv_quality"].isna().sum())
    if miss_cv:
        cnt.add("vOTUs with NO CheckV record", miss_cv)
    df["checkv_quality"] = df["checkv_quality"].fillna("Not-determined")
    df["CheckV_Provirus"] = df["CheckV_Provirus"].fillna("No")
    df["contig_length"] = pd.to_numeric(df["contig_length"],
                                        errors="coerce").fillna(0).astype(int)
    n_badlen = int((df["contig_length"] == 0).sum())
    if n_badlen:
        cnt.add("vOTUs with no usable CheckV length", n_badlen)

    # ---------------------------------------------------------------- 5
    print("[5] length flags (flags, not filters) ...")
    df["is_circular_checkv"] = (
        df["checkv_quality"].str.lower() == "complete").map({True: "Yes",
                                                            False: "No"})
    df["len_ge_5kb"] = (df["contig_length"] >= 5000).map({True: "Yes",
                                                          False: "No"})
    df["len_ge_10kb"] = (df["contig_length"] >= 10000).map({True: "Yes",
                                                            False: "No"})

    # ---------------------------------------------------------------- 6
    df["Base_Lifestyle"] = "Unknown"
    m_pro = (df["CheckV_Provirus"].str.strip() == "Yes") | \
            (df["Cluster_has_Provirus"] == "Yes")
    df.loc[m_pro, "Base_Lifestyle"] = "Provirus (CheckV/cluster)"
    print(f"[6] provirus evidence at this stage: {int(m_pro.sum())} vOTUs")

    cols = ["vOTU_ID", "contig_length", "checkv_quality", "completeness",
            "Zone", "Fraction", "Fraction_group", "MDA", "Mitomycin",
            "Quantified_own",
            "Cluster_size", "Cluster_has_Provirus", "CheckV_Provirus",
            "MAG_n", "MAG_ids",
            "is_circular_checkv", "len_ge_5kb", "len_ge_10kb",
            "Base_Lifestyle"]
    df = df[cols]
    df.to_csv(F_VOTU_01, sep="\t", index=False)

    # ---------------------------------------------------------------- report
    print("\n" + "=" * 70)
    print(" VERIFICATION")
    print("=" * 70)
    print(f" vOTUs in  : {n0}")
    print(f" vOTUs out : {len(df)}   (must be identical — this step never "
          f"adds or removes rows)")
    if len(df) != n0:
        print(" [FATAL] row count changed. A merge duplicated rows.")
        sys.exit(1)
    print("\n by Fraction_group:")
    for k, v in df["Fraction_group"].value_counts().items():
        print(f"   {v:>8}  {k}")
    print("\n length:")
    for c in ("len_ge_5kb", "len_ge_10kb", "is_circular_checkv"):
        print(f"   {int((df[c]=='Yes').sum()):>8}  {c}")
    cnt.report(total=n0)
    print("\n" + "=" * 70)
    print(f" Output: {F_VOTU_01}")


if __name__ == "__main__":
    main()
