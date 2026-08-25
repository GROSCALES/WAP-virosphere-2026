#!/usr/bin/env python3
"""
votus/05_add_ct3.py — Cenote-Taker3 hallmarks, topology and taxonomy.

Extracts hallmark genes (virion, replication, RDRP) and DTR/ITR topology 
from CT3 output. These represent two of the five independent lines of 
evidence used by the catalogue filter in step 06.

Consolidates taxonomy between CT3 and PhaBOX: the deepest assignment wins, 
and the winning source is recorded in Taxonomy_Source for full traceability.

Output: votu_05_ct3.tsv
"""

import os
import sys

import pandas as pd

from config import F_VOTU_04, F_VOTU_05, F_CT3, require, Counter_

NULLS = {"-", "", "nan", "none", "unknown", "unclassified", "na",
         "-_viruses;-_unclassified virus", "superkingdom:viruses"}

RANK_MAP = {"superkingdom": "-_", "clade": "-_", "kingdom": "k_",
            "phylum": "p_", "class": "c_", "order": "o_", "family": "f_",
            "subfamily": "sf_", "genus": "g_", "species": "s_"}

cnt = Counter_("COUNTERS — nothing is dropped in silence")


def is_null(x):
    return str(x).strip().lower() in NULLS


def phabox_to_ct3_style(tax):
    """PhaBOX writes 'rank:name;rank:name'; CT3 writes 'p_Name;c_Name'."""
    out = []
    for part in str(tax).split(";"):
        if ":" in part:
            rank, name = part.split(":", 1)
            out.append(f"{RANK_MAP.get(rank.strip().lower(), '-_')}{name.strip()}")
        else:
            out.append(f"-_{part.strip()}")
    return ";".join(out)


def depth(tax):
    return 0 if is_null(tax) else len([p for p in str(tax).split(";") if p.strip()])


def main():
    require(F_VOTU_04, F_CT3)
    print("=" * 70)
    print(" vOTU LEVEL 5 — Cenote-Taker3 + consolidated taxonomy")
    print("=" * 70)

    df = pd.read_csv(F_VOTU_04, sep="\t", dtype=str)
    n0 = len(df)
    print(f"[1] table in: {n0} vOTUs")

    # ---------------------------------------------------------------- 2
    print(f"[2] reading {os.path.basename(F_CT3)} ...")
    ct = pd.read_csv(F_CT3, sep="\t", dtype=str)
    print(f"    rows: {len(ct)}")

    want = {
        "input_name": "vOTU_ID",
        "end_feature": "CT3_Topology",
        "taxonomy_hierarchy": "CT3_Taxonomy",
        "virion_hallmark_count": "CT3_Virion_Hallmarks",
        "rep_hallmark_count": "CT3_Rep_Hallmarks",
        "RDRP_hallmark_count": "CT3_RDRP_Hallmarks",
    }
    missing = [c for c in want if c not in ct.columns]
    if missing:
        raise SystemExit(f"[SCHEMA ERROR] CT3 lacks columns: {missing}\n"
                         f"  found: {list(ct.columns)[:15]}")
    ct = ct[list(want)].rename(columns=want)
    ct["vOTU_ID"] = ct["vOTU_ID"].astype(str).str.strip()
    df["vOTU_ID"] = df["vOTU_ID"].astype(str).str.strip()

    dup = int(ct["vOTU_ID"].duplicated().sum())
    if dup:
        cnt.add("duplicate vOTU_IDs in the CT3 output (kept first)", dup)
        ct = ct.drop_duplicates(subset=["vOTU_ID"], keep="first")

    # --- merge audit, before merging (kept from the original) -----------
    n_match = int(df["vOTU_ID"].isin(ct["vOTU_ID"]).sum())
    n_orphan = int((~ct["vOTU_ID"].isin(df["vOTU_ID"])).sum())
    print(f"[3] merge audit")
    print(f"    vOTUs matched in CT3     : {n_match} / {n0} "
          f"({100.0*n_match/n0:.1f}%)")
    print(f"    CT3 rows with no vOTU    : {n_orphan}")
    if n_orphan:
        cnt.add("CT3 predictions whose contig is not in the catalogue",
                n_orphan)

    df = df.merge(ct, on="vOTU_ID", how="left")
    if len(df) != n0:
        print("[FATAL] CT3 merge changed the row count")
        sys.exit(1)

    cnt.add("vOTUs with NO CT3 record", n0 - n_match)
    df["CT3_Topology"] = df["CT3_Topology"].fillna("None")
    df["CT3_Taxonomy"] = df["CT3_Taxonomy"].fillna("Unknown")
    for c in ("CT3_Virion_Hallmarks", "CT3_Rep_Hallmarks",
              "CT3_RDRP_Hallmarks"):
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)

    df["CT3_Any_Hallmark"] = (
        (df["CT3_Virion_Hallmarks"] > 0) |
        (df["CT3_Rep_Hallmarks"] > 0) |
        (df["CT3_RDRP_Hallmarks"] > 0)).map({True: "Yes", False: "No"})
    df["CT3_Circular"] = df["CT3_Topology"].str.upper().isin(
        ["DTR", "ITR"]).map({True: "Yes", False: "No"})

    # ---------------------------------------------------------------- 4
    print("[4] consolidating taxonomy (CT3 vs PhaBOX) ...")

    def consolidate(r):
        ct3 = "" if is_null(r["CT3_Taxonomy"]) else str(r["CT3_Taxonomy"]).strip()
        pb = "" if is_null(r["PhaBOX_Taxonomy"]) else str(r["PhaBOX_Taxonomy"]).strip()
        d_ct3, d_pb = depth(ct3), depth(pb)
        if not ct3 and not pb:
            return pd.Series(["Unknown_Virus", "None", 0])
        if d_ct3 >= d_pb and ct3:
            return pd.Series([ct3, "CT3", d_ct3])
        return pd.Series([phabox_to_ct3_style(pb), "PhaBOX", d_pb])

    df[["Taxonomy", "Taxonomy_Source", "Taxonomy_Depth"]] = \
        df.apply(consolidate, axis=1)
    for src, n in df["Taxonomy_Source"].value_counts().items():
        cnt.add(f"taxonomy taken from {src}", int(n))

    df.to_csv(F_VOTU_05, sep="\t", index=False)

    # ---------------------------------------------------------------- report
    print("\n" + "=" * 70)
    print(" VERIFICATION")
    print("=" * 70)
    print(f" vOTUs in / out : {n0} / {len(df)}")

    print("\n CT3 evidence:")
    print(f"   with any hallmark   : "
          f"{int((df.CT3_Any_Hallmark=='Yes').sum())}")
    print(f"   virion hallmark > 0 : "
          f"{int((df.CT3_Virion_Hallmarks>0).sum())}")
    print(f"   rep hallmark > 0    : {int((df.CT3_Rep_Hallmarks>0).sum())}")
    print(f"   RDRP hallmark > 0   : {int((df.CT3_RDRP_Hallmarks>0).sum())}")
    print(f"   DTR/ITR (circular)  : {int((df.CT3_Circular=='Yes').sum())}")

    print("\n topology:")
    for k, v in df["CT3_Topology"].value_counts().head(6).items():
        print(f"   {v:>8}  {k}")

    print("\n taxonomy source:")
    for k, v in df["Taxonomy_Source"].value_counts().items():
        print(f"   {v:>8}  {k}   ({100.0*v/n0:.1f}%)")

    print("\n top taxa (consolidated):")
    for k, v in df["Taxonomy"].value_counts().head(10).items():
        print(f"   {v:>7}  {str(k)[:85]}")

    print("\n independent evidence available per vOTU "
          "(what the step-06 filter will use):")
    ev = pd.DataFrame({
        "VIBRANT": df.Vibrant_predicted == "Yes",
        "PhaBOX": df.PhaBOX_Pred != "Unknown",
        "CT3_hallmark": df.CT3_Any_Hallmark == "Yes",
        "CheckV_quality": df.checkv_quality.str.lower().isin(
            ["medium-quality", "high-quality", "complete"]),
        "Provirus": df.CheckV_Provirus.str.strip() == "Yes",
    })
    for c in ev.columns:
        print(f"   {int(ev[c].sum()):>8}  {c}")
    n_ev = ev.sum(axis=1)
    print("\n number of independent lines of evidence:")
    for k in range(6):
        n = int((n_ev == k).sum())
        if n:
            print(f"   {n:>8}  {k} source(s)")

    cnt.report(total=n0)
    print("\n" + "=" * 70)
    print(f" Output: {F_VOTU_05}")


if __name__ == "__main__":
    main()
