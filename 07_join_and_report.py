#!/usr/bin/env python3
"""
07_join_and_report.py — Joins the vOTU and gene branches to produce final reporting tables.

Inputs:
  votus/06 -> Final vOTU catalog flags
  genes/05 -> Gene classification (Defence, AVG, Other)
  genes/06 -> Harmonised category map

Outputs:
  07_long_table_FINAL.tsv : One row per Defence/AVG gene with full vOTU metadata.
  07_votu_with_counts.tsv : The vOTU table annotated with n_orfs, n_AVG, n_Defence.

Design Principles:
1. Single Source of Truth for Counts: Per-vOTU gene counts are derived strictly 
   by aggregating the gene table here. They are never computed independently 
   to guarantee absolute synchrony between the long and wide formats.
2. System-Level Categories: Categories reflect the defence SYSTEM a gene belongs to, 
   not the isolated biochemical activity of the individual protein. For example, 
   an RM system comprises both a nucleolytic R subunit and a methyltransferase M subunit; 
   both are filed under the system's primary mechanism (Nucleic acid degradation).
3. Downstream Filtering: Filtering happens at reporting time. The long table carries 
   boolean flags (Passes_catalog, Conservative_10kb, passes_context) so sensitivity 
   analyses can be seamlessly executed downstream in R.
"""

import os
from collections import defaultdict, Counter

import pandas as pd

from config import (
    F_GENE_TABLE_CLASS, F_VOTU_06, F_LONG_TABLE, OUT_DIR,
    DEFENCE_CATEGORIES, require, Counter_,
)

F_VOTU_COUNTS = os.path.join(OUT_DIR, "07_votu_with_counts.tsv")
F_CAT_MAP = os.path.join(OUT_DIR, "06_category_map.tsv")

# Numbers reported in the current manuscript, for comparison only.
OLD_AMG_TOTAL = 1441
OLD_DEFENCE_TOTAL = 1964

cnt = Counter_("COUNTERS — nothing is dropped in silence")


def load_category_map(path):
    """raw defence label -> harmonised category."""
    m = {}
    if not os.path.exists(path):
        raise SystemExit(f"[MISSING] {path}\n"
                         f"  run genes/06_harmonise_categories.py first")
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i = {c: j for j, c in enumerate(header)}
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) > max(i["category_raw"], i["category"]):
                m[c[i["category_raw"]]] = c[i["category"]]
    return m


def main():
    require(F_GENE_TABLE_CLASS, F_VOTU_06)
    print("=" * 74)
    print(" JOIN — genes x vOTUs, categories harmonised, filtering left to R")
    print("=" * 74)

    catmap = load_category_map(F_CAT_MAP)
    print(f"[0] category map: {len(catmap)} raw labels")
    n_unmapped_labels = sum(1 for v in catmap.values() if v == "UNMAPPED")
    if n_unmapped_labels:
        cnt.add("raw labels still UNMAPPED in the category map",
                n_unmapped_labels)

    print("[1] vOTU catalogue ...")
    votu = pd.read_csv(F_VOTU_06, sep="\t", dtype=str, low_memory=False)
    n_votu, n_pass = len(votu), int((votu["Passes_catalog"] == "Yes").sum())
    print(f"    vOTUs: {n_votu}   in final catalogue: {n_pass}")

    print("[2] gene table (streaming to keep memory flat) ...")
    with open(F_GENE_TABLE_CLASS) as fh:
        gh = fh.readline().rstrip("\n").split("\t")
    gi = {c: i for i, c in enumerate(gh)}
    for c in ("protein_id", "votu_id", "gene_class", "category",
              "class_source", "passes_context", "aux_score"):
        if c not in gi:
            raise SystemExit(f"[SCHEMA ERROR] gene table lacks '{c}'")

    n_orfs, n_avg, n_def = defaultdict(int), defaultdict(int), defaultdict(int)
    keep, n_excluded = [], 0

    with open(F_GENE_TABLE_CLASS) as fh:
        next(fh)
        for line in fh:
            c = line.rstrip("\n").split("\t")
            v = c[gi["votu_id"]]
            n_orfs[v] += 1
            klass = c[gi["gene_class"]]
            if klass == "Other":
                continue

            raw = c[gi["category"]]
            if klass == "Defence":
                final = catmap.get(raw, "UNMAPPED")
                if final.startswith("EXCLUDED"):
                    c[gi["gene_class"]] = "Excluded_not_defence"
                    n_excluded += 1
                    cnt.add(f"reclassified out of Defence: {raw}")
                else:
                    n_def[v] += 1
                    if final == "UNMAPPED":
                        cnt.add(f"defence label with no category: {raw}")
            else:                       # AVG keeps its KEGG category
                final = raw or "Unclassified"
                n_avg[v] += 1

            c.append(final)
            keep.append(c)

    print(f"    ORFs read        : {sum(n_orfs.values())}")
    print(f"    Defence + AVG    : {len(keep)}")
    print(f"    of which excluded: {n_excluded}")

    genes = pd.DataFrame(keep, columns=gh + ["category_final"])

    print("[3] deriving per-vOTU counts from the gene table ...")
    votu["n_orfs"] = votu["vOTU_ID"].map(n_orfs).fillna(0).astype(int)
    votu["n_AVG"] = votu["vOTU_ID"].map(n_avg).fillna(0).astype(int)
    votu["n_Defence"] = votu["vOTU_ID"].map(n_def).fillna(0).astype(int)
    votu.to_csv(F_VOTU_COUNTS, sep="\t", index=False)

    print("[4] joining vOTU metadata onto the genes ...")
    meta_cols = [c for c in [
        "vOTU_ID", "contig_length", "checkv_quality", "completeness",
        "Zone", "Fraction", "Fraction_group", "MDA", "Mitomycin",
        "Quantified_own", "MAG_n", "MAG_ids", "MAG_taxonomy",
        "Taxonomy", "Taxonomy_Source", "Host_Prediction", "Host_Source",
        "Lifestyle", "Lifestyle_Source", "N_evidence", "Evidence_sources",
        "Passes_catalog", "Conservative_10kb", "len_ge_5kb", "len_ge_10kb",
        "TPM_total", "n_samples_detected", "Detected_quantitative",
        "n_orfs", "n_AVG", "n_Defence",
    ] if c in votu.columns]
    tpm_cols = [c for c in votu.columns
                if c.startswith("TPM_") and c not in meta_cols]

    long = genes.merge(votu[meta_cols + tpm_cols],
                       left_on="votu_id", right_on="vOTU_ID", how="left")
    miss = int(long["Passes_catalog"].isna().sum())
    if miss:
        cnt.add("gene rows with no vOTU metadata", miss)
    long["Passes_catalog"] = long["Passes_catalog"].fillna("No")
    long.to_csv(F_LONG_TABLE, sep="\t", index=False)

    # ---------------------------------------------------------------- report
    fin = long[(long["Passes_catalog"] == "Yes") &
               (long["gene_class"].isin(["AVG", "Defence"]))]
    ctx = fin[fin["passes_context"] == "1"]
    cons = fin[fin["Conservative_10kb"] == "Yes"]

    def n(d, k):
        return int((d["gene_class"] == k).sum())

    print("\n" + "=" * 74)
    print(" RESULT — final catalogue only")
    print("=" * 74)
    print(f"{'subset':<34}{'AVG':>10}{'Defence':>10}{'vOTUs':>9}")
    print("-" * 74)
    for label, d in [("all genes in catalogue", fin),
                     ("+ context aux_score <= 3", ctx),
                     ("+ conservative >= 10 kb", cons)]:
        print(f"{label:<34}{n(d,'AVG'):>10}{n(d,'Defence'):>10}"
              f"{d['votu_id'].nunique():>9}")

    print("\n" + "-" * 74)
    print(" COMPARISON WITH THE CURRENT MANUSCRIPT")
    print("-" * 74)
    for label, old, new in [("AVG / AMG", OLD_AMG_TOTAL, n(fin, "AVG")),
                            ("Defence", OLD_DEFENCE_TOTAL, n(fin, "Defence"))]:
        d = new - old
        print(f"   {label:<12} manuscript {old:>6}   new {new:>6}   "
              f"{d:+6} ({100.0*d/old:+.1f}%)")

    print("\n" + "-" * 74)
    print(" DEFENCE BY HARMONISED CATEGORY (final catalogue)")
    print(" categories are SYSTEM-level: an RM system's methyltransferase")
    print(" subunit is filed under the system, not under 'modification'")
    print("-" * 74)
    d = fin[fin["gene_class"] == "Defence"]
    tot = len(d)
    for c in DEFENCE_CATEGORIES:
        k = int((d["category_final"] == c).sum())
        if k:
            print(f"   {k:>6}  {c:<38} ({100.0*k/tot:.1f}%)")
    other = int((~d["category_final"].isin(DEFENCE_CATEGORIES)).sum())
    if other:
        print(f"   {other:>6}  (unmapped / other)")

    print("\n AVG BY KEGG CATEGORY (final catalogue):")
    a = fin[fin["gene_class"] == "AVG"]
    for c, k in a["category_final"].value_counts().head(12).items():
        print(f"   {k:>6}  {c}")

    print("\n" + "-" * 74)
    print(" BUTTERFLY — % of vOTUs encoding >=1 gene of each class")
    print("-" * 74)
    vf = votu[votu["Passes_catalog"] == "Yes"]
    rows = []
    for frac in ["CB22", "VIRFC", "GLI", "LIS"]:
        s = vf[vf["Fraction"] == frac]
        if len(s):
            rows.append({"Fraction": frac, "vOTUs": len(s),
                         "%with_AVG": round(100.0*(s["n_AVG"] > 0).mean(), 1),
                         "%with_Defence": round(100.0*(s["n_Defence"] > 0).mean(), 1)})
    print(pd.DataFrame(rows).to_string(index=False))
    print("\n   GLI vs CB22 is the clean contrast: both non-MDA, neither")
    print("   mitomycin-induced. Each row is a 2x2 Fisher test.")

    print("\n categories x culture group (the Figure 3 matrix):")
    tab = pd.crosstab(d["category_final"], d["Fraction_group"])
    print(tab.to_string())

    cnt.report()
    print("\n" + "=" * 74)
    print(f" Long table : {F_LONG_TABLE}")
    print(f" vOTU table : {F_VOTU_COUNTS}")
    print(" R plots category_final; filter on Passes_catalog,")
    print(" Conservative_10kb and passes_context there, not here.")


if __name__ == "__main__":
    main()
