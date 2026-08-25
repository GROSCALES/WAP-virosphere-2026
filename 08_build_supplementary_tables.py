#!/usr/bin/env python3
# =============================================================================
# 08_build_supplementary_tables.py
#
# Builds the two clean, publication-ready supplementary tables from the
# pipeline's own source of truth (07_votu_with_counts.tsv,
# 07_long_table_FINAL.tsv). This is deliberately the LAST step of the
# Python pipeline, not something built in R: R is for stats and figures,
# and building "what counts as the final table" in a second language on a
# second machine is exactly the kind of duplication that produced the
# three-different-AMG-totals bug in the old Table S3. One producer, one
# source of truth.
#
#   Table_S2_vOTU_catalogue.tsv     one row per vOTU that passes the catalog
#   Table_S3_AVG_Defence_genes.tsv  one row per AVG/Defence gene that
#                                   contributes to those vOTUs' n_AVG/n_Defence
#
# WHAT IS DELIBERATELY LEFT OUT (and why)
#   - Per-sample TPM_ZONE_FRACTION columns in Table S2: collapsed into a
#     single CPM_total. The full abundance matrix is useful but belongs in
#     its own file if you want to publish it (say so and I'll add a
#     Table_S_abundance_matrix.tsv), not bloating the main catalogue table.
#   - Vetoed / excluded genes in Table S3: the veto logic and what got
#     excluded and why lives in the pipeline code (05_classify.py,
#     config.py's CORE_FUNCTION_REJECT), which is going to GitHub/Zenodo.
#     Duplicating that trace into a giant supplementary table would make
#     S3 unreadable for the one thing a reader actually needs it for:
#     "what genes make up the numbers in Figure 3".
#
# CHECKS, so a wrong column pick fails loudly instead of silently
#   - row count of Table_S3 grouped by (vOTU_ID, gene_class) sums back to
#     exactly n_AVG / n_Defence in Table_S2 for every vOTU (the same
#     invariant 99_qc.py already checks upstream; re-checked here because
#     this script re-derives its own counts independently as a second,
#     cheap verification).
#   - no duplicate protein_id.
# =============================================================================

import os
import pandas as pd
from pathlib import Path

OUT_DIR = Path(os.environ.get(
    "WAP_BASE", Path.home() / "gabi/vOTUs_total")) / "pipeline_v2" / "out"

votu_path = OUT_DIR / "07_votu_with_counts.tsv"
long_path = OUT_DIR / "07_long_table_FINAL.tsv"

votu = pd.read_csv(votu_path, sep="\t", low_memory=False)
long = pd.read_csv(long_path, sep="\t", low_memory=False)

def require(df, cols, label):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise SystemExit(f"{label}: missing columns {missing} -- "
                         f"header must have changed, check the table.")

# =============================================================================
# TABLE S2 — vOTU catalogue
# =============================================================================
S2_TPM_COLS = ["TPM_DE_CB22", "TPM_DE_GLI", "TPM_EL_CB22", "TPM_EL_GLI",
              "TPM_MB_CB22", "TPM_MB_GLI", "TPM_OR_CB22", "TPM_OR_GLI"]

s2_source_cols = [
    "vOTU_ID", "Fraction", "Fraction_group", "Zone", "contig_length",
    "completeness", "checkv_quality", "N_evidence", "Lifestyle",
    "Lifestyle_Source", "Taxonomy", "Taxonomy_Source", "Host_Prediction",
    "Host_Source", "n_orfs", "n_AVG", "n_Defence", "Conservative_10kb",
    "MDA", "Mitomycin", "Quantified_own",
] + S2_TPM_COLS

require(votu, ["Passes_catalog"] + s2_source_cols, "07_votu_with_counts.tsv")

s2 = votu.loc[votu["Passes_catalog"] == "Yes", s2_source_cols].copy()
s2["CPM_total"] = s2[S2_TPM_COLS].sum(axis=1)
s2 = s2.drop(columns=S2_TPM_COLS)

s2 = s2.rename(columns={
    "vOTU_ID": "vOTU_ID", "Fraction": "Fraction", "Fraction_group": "Origin",
    "Zone": "Zone", "contig_length": "Contig_length_bp",
    "completeness": "Completeness_pct", "checkv_quality": "CheckV_quality",
    "N_evidence": "N_evidence_lines", "Lifestyle": "Lifestyle",
    "Lifestyle_Source": "Lifestyle_source", "Taxonomy": "Viral_taxonomy",
    "Taxonomy_Source": "Taxonomy_source", "Host_Prediction": "Predicted_host",
    "Host_Source": "Host_prediction_source", "n_orfs": "N_ORFs",
    "n_AVG": "N_AVG_genes", "n_Defence": "N_Defence_genes",
    "Conservative_10kb": "Meets_10kb_or_circular",
    "MDA": "MDA_amplified", "Mitomycin": "Mitomycin_induced",
    "Quantified_own": "Own_reads_quantified", "CPM_total": "CPM_total",
})

# =============================================================================
# TABLE S3 — AVG + Defence gene catalogue
# =============================================================================
# df_in_complete_system / pl_in_complete_system were computed back in
# 02_add_evidence.py and then never used anywhere. They are the only columns
# that distinguish "a component of a defence system that was actually
# detected in full" from "an isolated domain hit", and the annotation audit
# turned that distinction into a finding worth reporting: only ~6% of the
# defence genes here belong to a complete system, and for the CRISPR-Cas
# category it was 0 of 176. That is the expected result on phage contigs,
# which carry single genes and fragments rather than whole bacterial defence
# operons -- but a reviewer will ask, so the evidence ships with the table.
s3_source_cols = [
    "protein_id", "vOTU_ID", "gene_class", "category_final",
    "ko", "ko_name", "pfam", "pfam_name", "df_system", "pl_system",
    "evidence_sources", "aux_score", "passes_context", "specialist_vetoed",
    "df_in_complete_system", "pl_in_complete_system",
]
require(long, ["Passes_catalog", "gene_class"] + s3_source_cols,
       "07_long_table_FINAL.tsv")

s3 = long.loc[
    (long["Passes_catalog"] == "Yes") & (long["gene_class"].isin(["AVG", "Defence"])),
    s3_source_cols
].copy()

s3 = s3.rename(columns={
    "protein_id": "Gene_ID", "vOTU_ID": "vOTU_ID", "gene_class": "Gene_class",
    "category_final": "Category", "ko": "KO", "ko_name": "KO_name",
    "pfam": "Pfam", "pfam_name": "Pfam_name",
    "df_system": "DefenseFinder_raw_system", "pl_system": "PADLOC_raw_system",
    "evidence_sources": "Evidence_sources", "aux_score": "DRAMv_aux_score",
    "passes_context": "Aux_score_le_3", "specialist_vetoed": "Core_function_vetoed",
    "df_in_complete_system": "In_complete_DefenseFinder_system",
    "pl_in_complete_system": "In_complete_PADLOC_system",
})
s3["In_complete_system"] = (
    (s3["In_complete_DefenseFinder_system"].astype(str) == "1") |
    (s3["In_complete_PADLOC_system"].astype(str) == "1")
).map({True: "Yes", False: "No"})

# =============================================================================
# CHECKS
# =============================================================================
dup = s3["Gene_ID"].duplicated().sum()
if dup:
    raise SystemExit(f"Table S3: {dup} duplicate Gene_ID -- stop, do not export.")

recount = (s3.groupby(["vOTU_ID", "Gene_class"]).size()
          .unstack(fill_value=0)
          .reindex(columns=["AVG", "Defence"], fill_value=0))
check = s2.set_index("vOTU_ID")[["N_AVG_genes", "N_Defence_genes"]].join(recount, how="left").fillna(0)
mismatch = check[(check["N_AVG_genes"] != check.get("AVG", 0)) |
                (check["N_Defence_genes"] != check.get("Defence", 0))]
if len(mismatch):
    print(f"WARNING: {len(mismatch)} vOTUs where Table S2's N_AVG/N_Defence "
         f"don't match Table S3's own row counts -- inspect before shipping:")
    print(mismatch.head(10))
else:
    print(f"OK: Table S3 row counts reconcile exactly with Table S2's "
         f"N_AVG_genes/N_Defence_genes for all {len(check)} vOTUs.")

# =============================================================================
# EXPORT
# =============================================================================
s2_out = OUT_DIR / "Table_S2_vOTU_catalogue.tsv"

# ==========================================================================
# PRESENTATION LABELS — applied here and nowhere else
# ==========================================================================
# The manuscript uses the phage-biology convention: a PHAGE is virulent or
# temperate, a CYCLE is lytic or lysogenic. The pipeline's internal value is
# "Lytic" because that is what VIBRANT and PhaBOX return, and it is kept
# internally so the code matches the tools it wraps.
#
# The manuscript also says AMG (auxiliary metabolic gene) where the code says
# AVG (auxiliary viral gene) -- the same set of genes, named differently
# because that is what the variable was called while the pipeline was written.
#
# Both are relabelled HERE, at the boundary between the pipeline and the
# reader, and AFTER the reconciliation check above, which needs the internal
# names. Nothing downstream reads these two files except people.
s2["Lifestyle"] = s2["Lifestyle"].replace({"Lytic": "Virulent"})
s2 = s2.rename(columns={"N_AVG_genes": "N_AMG_genes"})
s3["Gene_class"] = s3["Gene_class"].replace({"AVG": "AMG"})

s3_out = OUT_DIR / "Table_S3_AVG_Defence_genes.tsv"
s2.to_csv(s2_out, sep="\t", index=False)
s3.to_csv(s3_out, sep="\t", index=False)

print(f"\nTable S2: {len(s2)} vOTUs, {len(s2.columns)} columns -> {s2_out}")
print(f"Table S3: {len(s3)} genes ({(s3.Gene_class=='AMG').sum()} AMG, "
     f"{(s3.Gene_class=='Defence').sum()} Defence), {len(s3.columns)} columns -> {s3_out}")
_d = s3[s3.Gene_class == "Defence"]
_c = (_d.In_complete_system == "Yes").sum()
print(f"\n  of the {len(_d)} defence genes, {_c} ({100.0*_c/len(_d):.1f}%) belong to a "
     f"COMPLETE system.\n  Low is expected on phage contigs; report the number rather "
     f"than letting a reviewer find it.")
print("  by category:")
for cat, g in _d.groupby("Category"):
    k = (g.In_complete_system == "Yes").sum()
    print(f"    {len(g):>5}  {k:>4} complete ({100.0*k/len(g):5.1f}%)  {cat}")
print("\nTable S2 columns:", list(s2.columns))
print("Table S3 columns:", list(s3.columns))
