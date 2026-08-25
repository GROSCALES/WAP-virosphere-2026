#!/usr/bin/env python3
"""
99_qc.py — Quality Control gate and manuscript statistics generator.

This script serves two primary functions:

  1. STRUCTURAL CONSISTENCY: Enforces data integrity across the entire pipeline.
     Asserts that row counts are strictly preserved across merges, IDs are 
     unique, per-vOTU derived counts perfectly match the explicit aggregation 
     of the gene table, length/evidence flags are logically consistent, and 
     the category map is completely resolved. Exits non-zero on any failure, 
     acting as a strict gate for automated execution (e.g., Makefiles).

  2. MANUSCRIPT METRICS: Compiles all summary statistics required for the 
     manuscript into a single, centralized output block. This guarantees that 
     all reported numbers are deterministically generated directly from the 
     final tables, preventing manual tallying errors.

Note: This is a read-only execution; it modifies no pipeline outputs.
"""

import os
import sys

import pandas as pd

from config import (
    OUT_DIR, F_VOTU_01, F_VOTU_06, F_GENE_TABLE_CLASS, F_LONG_TABLE,
    F_QC_REPORT, DEFENCE_CATEGORIES, LEN_MAIN_MIN, LEN_CONSERVATIVE_MIN,
    EVALUE_MAX, CTX_WINDOW_ORFS, require,
)

F_VOTU_COUNTS = os.path.join(OUT_DIR, "07_votu_with_counts.tsv")
F_CAT_MAP = os.path.join(OUT_DIR, "06_category_map.tsv")

_fail = []
_out = []


def say(s=""):
    print(s)
    _out.append(s)


def check(name, ok, detail=""):
    say(f"   [{'PASS' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        _fail.append(name)


def main():
    require(F_VOTU_01, F_VOTU_06, F_GENE_TABLE_CLASS, F_LONG_TABLE,
            F_VOTU_COUNTS, F_CAT_MAP)

    say("=" * 74)
    say(" QC — consistency checks")
    say("=" * 74)

    v01 = pd.read_csv(F_VOTU_01, sep="\t", usecols=["vOTU_ID"], dtype=str)
    votu = pd.read_csv(F_VOTU_COUNTS, sep="\t", dtype=str, low_memory=False)
    long = pd.read_csv(F_LONG_TABLE, sep="\t", dtype=str, low_memory=False)
    cmap = pd.read_csv(F_CAT_MAP, sep="\t", dtype=str)

    for c in ("n_orfs", "n_AVG", "n_Defence", "contig_length"):
        if c in votu.columns:
            votu[c] = pd.to_numeric(votu[c], errors="coerce").fillna(0)

    # ---- 1 vOTU rows preserved end to end -------------------------------
    say("\n vOTU branch")
    check("row count preserved from step 01 to step 06",
          len(v01) == len(votu), f"{len(v01)} -> {len(votu)}")
    check("no duplicated vOTU_ID",
          not votu["vOTU_ID"].duplicated().any())

    # ---- 2 gene table ----------------------------------------------------
    say("\n gene branch")
    n_gene_rows = sum(1 for _ in open(F_GENE_TABLE_CLASS)) - 1
    check("no duplicated protein_id in the long table",
          not long["protein_id"].duplicated().any(),
          f"{long['protein_id'].nunique()} unique of {len(long)}")
    check("every long-table gene has vOTU metadata",
          long["Passes_catalog"].notna().all())

    # ---- 3 derived counts really are derived -----------------------------
    say("\n derived counts (the bug that forced the old manual patches)")
    agg = (long[long["gene_class"].isin(["AVG", "Defence"])]
           .groupby(["votu_id", "gene_class"]).size().unstack(fill_value=0))
    for cls, col in (("AVG", "n_AVG"), ("Defence", "n_Defence")):
        want = votu.set_index("vOTU_ID")[col]
        got = agg[cls] if cls in agg.columns else pd.Series(dtype=int)
        joined = want.to_frame("table").join(got.rename("long"), how="left")
        joined["long"] = joined["long"].fillna(0)
        bad = int((joined["table"] != joined["long"]).sum())
        check(f"{col} equals the long-table aggregation", bad == 0,
              f"{bad} mismatching vOTUs")

    # ---- 4 flags nested correctly ---------------------------------------
    say("\n flags")
    pc = votu["Passes_catalog"] == "Yes"
    c10 = votu["Conservative_10kb"] == "Yes"
    check("Conservative_10kb is a subset of Passes_catalog",
          int((c10 & ~pc).sum()) == 0)
    check("every catalogue vOTU has length >= threshold, is circular, "
          "or was rescued",
          int((pc & (votu["contig_length"] < LEN_MAIN_MIN) &
               (votu.get("Rescued_short", "No") != "Yes") &
               (votu["is_circular_checkv"] != "Yes") &
               (votu["CT3_Circular"] != "Yes")).sum()) == 0)

    # ---- 5 category map --------------------------------------------------
    say("\n categories")
    unmapped = int((cmap["category"] == "UNMAPPED").sum())
    check("category map has no UNMAPPED label", unmapped == 0,
          f"{unmapped} labels")
    d = long[(long["Passes_catalog"] == "Yes") &
             (long["gene_class"] == "Defence")]
    bad_cat = int((~d["category_final"].isin(DEFENCE_CATEGORIES)).sum())
    check("every defence gene has a category from the declared vocabulary",
          bad_cat == 0, f"{bad_cat} genes")

    # ==================================================================
    say("\n" + "=" * 74)
    say(" NUMBERS FOR THE MANUSCRIPT")
    say(" copy from here; do not re-derive by hand")
    say("=" * 74)

    fin = votu[pc]
    lf = long[long["Passes_catalog"] == "Yes"]
    la = lf[lf["gene_class"] == "AVG"]
    ld = lf[lf["gene_class"] == "Defence"]

    say("\n CATALOGUE")
    say(f"   vOTUs assembled                  {len(votu):>8}")
    say(f"   vOTUs in the final catalogue     {len(fin):>8}")
    say(f"   conservative subset (>= {LEN_CONSERVATIVE_MIN//1000} kb)     "
        f"{int(c10.sum()):>8}")
    ne = pd.to_numeric(fin["N_evidence"], errors="coerce")
    say(f"   supported by >= 2 evidence lines {int((ne >= 2).sum()):>8}  "
        f"({100.0*(ne >= 2).mean():.1f}%)")

    say("\n LIFESTYLE (final catalogue)")
    for k, n in fin["Lifestyle"].value_counts().items():
        say(f"   {n:>7}  {k}   ({100.0*n/len(fin):.1f}%)")

    say("\n LIFESTYLE BY FRACTION  (% temperate; GLI vs CB22 is the")
    say(" clean contrast: both non-MDA, neither mitomycin-induced)")
    for frac in ["CB22", "VIRFC", "GLI", "LIS"]:
        s = fin[fin["Fraction"] == frac]
        if len(s):
            t = int(s["Lifestyle"].str.startswith("Temperate").sum())
            say(f"   {frac:<6} {t:>5} / {len(s):<6} ({100.0*t/len(s):5.1f}%)")

    say("\n GENES (final catalogue)")
    say(f"   auxiliary viral genes (AVG)      {len(la):>8}")
    say(f"   defence-associated genes         {len(ld):>8}")
    say(f"   vOTUs encoding >= 1 AVG          {int((fin['n_AVG'] > 0).sum()):>8}")
    say(f"   vOTUs encoding >= 1 defence gene {int((fin['n_Defence'] > 0).sum()):>8}")

    say("\n DEFENCE BY CATEGORY (system-level; see the figure legend note)")
    for c in DEFENCE_CATEGORIES:
        n = int((ld["category_final"] == c).sum())
        if n:
            say(f"   {n:>7}  {c:<38} ({100.0*n/len(ld):.1f}%)")

    say("\n PER-vOTU DENSITY (raw counts are confounded by group size:")
    say(" the culture-independent group has ~10x more vOTUs)")
    for g in ["Culture-independent", "Culture-dependent"]:
        s = fin[fin["Fraction_group"] == g]
        if len(s):
            say(f"   {g:<22} {len(s):>6} vOTUs | "
                f"{s['n_Defence'].sum()/len(s):.3f} defence genes/vOTU | "
                f"{s['n_AVG'].sum()/len(s):.3f} AVG/vOTU")

    say("\n SENSITIVITY (report both; the conclusion should not change)")
    for label, sub in [("all genes", lf),
                       (f"aux_score <= 3", lf[lf["passes_context"] == "1"]),
                       (f">= {LEN_CONSERVATIVE_MIN//1000} kb or circular",
                        lf[lf["Conservative_10kb"] == "Yes"])]:
        say(f"   {label:<26} AVG {int((sub.gene_class=='AVG').sum()):>6}   "
            f"Defence {int((sub.gene_class=='Defence').sum()):>6}")

    say("\n PARAMETERS TO STATE IN METHODS")
    say(f"   HMM e-value cutoff               {EVALUE_MAX}  (as VIBRANT)")
    say(f"   context window                   +/- {CTX_WINDOW_ORFS} ORFs")
    say(f"   auxiliary score retained         <= 3  (DRAM-v default)")
    say(f"   minimum contig length            {LEN_MAIN_MIN} bp")
    say(f"   conservative subset              {LEN_CONSERVATIVE_MIN} bp or circular")

    with open(F_QC_REPORT, "w") as fh:
        fh.write("\n".join(_out) + "\n")

    say("\n" + "=" * 74)
    if _fail:
        say(f" {len(_fail)} CHECK(S) FAILED — do not plot yet:")
        for f in _fail:
            say(f"   - {f}")
        say("=" * 74)
        sys.exit(1)
    say(" ALL CHECKS PASSED")
    say(f" Report written to {F_QC_REPORT}")
    say("=" * 74)


if __name__ == "__main__":
    main()
