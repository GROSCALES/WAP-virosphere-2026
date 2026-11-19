# WAP-virosphere

Code and curated files behind **"Viral compartment, not viral lifestyle, predicts
phage-encoded defence strategies in the Southern Ocean"**.

This is the analysis record for one paper, not a maintained tool. Raw inputs are not
included, so the pipeline does not run end to end from here; what it gives is a complete
account of how every number and every figure was obtained. Read it alongside the
Materials and Methods and the Supplementary Methods, whose section numbers are cited
throughout.

Gabriel Roscales · Water Microbiology Related to Health group (MARS), University of
Barcelona. Archived at submission, November 2026.

## The paper

Polar viromics has focused almost exclusively on free viral particles, leaving the
integrated prophage reservoir largely unexplored. We combined four fractions — free
viral particles (VIR), the cellular fraction (BAC), bacterial cultures (CUL) and their
mitomycin C-induced lysates (LIS) — with hybrid short- and long-read assembly across
four stations, yielding 11,100 non-redundant prokaryotic vOTUs. Taxonomy, lifestyle and
functional content all partitioned by culture dependence, and the two classes of
accessory gene went in opposite directions: auxiliary metabolic genes are denser in the
culture-independent fractions, complete defence loci in the culture-dependent ones.
Defence loci resolve along two axes, self-protection and exclusion, and the balance
between them follows the fraction a virus was recovered from rather than the phage's
lifestyle.

### Key numbers

From `99_qc.py`. The manuscript copies these rather than re-deriving them.

| | |
|---|---|
| representative sequences screened | 59,469 |
| met the length criterion | 11,671 |
| eukaryotic viruses and virophages removed | 571 |
| vOTUs in the final catalogue | 11,100 |
| supported by ≥ 2 lines of evidence besides VirSorter2 | 10,671 (96.1%) |
| supported by exactly one | 429 (3.9%) |
| medium CheckV quality or better | 937 (8.4%) |
| conservative subset (≥ 10 kb or complete) | 3,810 |
| near-complete subset (gene-sharing network) | 363 |
| auxiliary metabolic genes (AMG) | 1,441 |
| defence-associated genes | 1,964 |
| reconstructed defence loci | 1,718 |
| complete defence loci | 808 |
| complete loci on the two antiviral axes | 773 |

The evidence requirement is descriptive and removes no vOTU: membership in the catalogue
is decided by length, topology and domain.

Defence genes by category, gene level:

| Category | n | % |
|---|---|---|
| Nucleic acid degradation | 840 | 42.8 |
| Other / uncharacterised systems | 386 | 19.7 |
| Abortive infection | 383 | 19.5 |
| Genome modification | 318 | 16.2 |
| Host manipulation (phage-encoded) | 31 | 1.6 |
| CRISPR-Cas | 6 | 0.3 |

Complete defence loci by axis:

| Axis | n | % |
|---|---|---|
| Self-protection | 599 | 77.5 |
| Exclusion | 174 | 22.5 |

Gene-level and locus-level distributions differ because a locus whose
restriction–modification call is rescinded for lack of an endonuclease is recounted as an
orphan methyltransferase under *Genome modification*.

---

## How it is organised

Two parallel branches, joined at the end.

```
config.py                          paths, thresholds, search patterns,
                                   category vocabulary
votus/                             one row per vOTU
  01_build_base.py                 04_add_phabox_gtdb.py
  02_add_vibrant.py                05_add_ct3.py
  03_add_coverm.py                 06_filter_catalog.py
genes/                             one row per ORF
  01_build_gene_table.py           04_score_context.py
  02_add_evidence.py               05_classify.py
  03_kegg_categories.py            06_harmonise_categories.py
07_join_and_report.py              joins the two branches
reconstruir_loci.py                genes -> defence loci
08_build_supplementary_tables.py   the vOTU and gene tables of the paper
08_defensas_votus.sh               PADLOC + DefenseFinder screening
98_audit_annotations.py            every gene against its own annotation
99_qc.py                           consistency checks + manuscript numbers
data/
  mapa_v5.tsv                      curated family knowledge base
  mapa_v5_EN.tsv                   English translation, for reading only
R/
  01_stats.R  02_Fig_2.R           community statistics and Figure 2
  03_stats_loci.R  04_stats_fig3.R statistics behind Figure 3
  05_Fig_3.R                       Figure 3
  06_Fig_4.R  07_Fig_5.R  08_Fig_6.R   genomic maps
  09_supplementary_statistics.R    Tables S3 and S4 from the published data
  palette_paper.R                  colours and vocabulary, Figures 2 and 3
  palette_maps.R                   colours and vocabulary, Figures 4, 5 and 6
requirements.txt  environment.txt  README.md  LICENSE
```

All paths, thresholds and category definitions live in `config.py`. If any of them
appears elsewhere, that is a bug.

### Nomenclature: the paper's words and the code's words

Four terms were renamed while the manuscript was written, after the pipeline had already
produced the results. The internal names were left alone, so that the code here is the
code the results came from.

| In the paper | In the code | Where the relabel happens |
|---|---|---|
| **AMG** (auxiliary metabolic gene) | `AVG`: `gene_class == "AVG"`, `n_AVG`, `Table_S3_AVG_Defence_genes.tsv` | `08_build_supplementary_tables.py` |
| **Virulent** | `Lytic`, which is what VIBRANT and PhaBOX return | `08_build_supplementary_tables.py` |
| **Self-protection** / **Exclusion** | `1 Contra-defensa` / `2 Defensa` | `reconstruir_loci.py`, via `data/mapa_v5.tsv` |
| **culture dependence** | `Origin`, `Fraction_group` | not relabelled |

`Origin` is the only one not renamed at the boundary. It was the working name for the
culture-dependent / culture-independent split, and it survives as a column in Data S1 and
as the string the figure scripts match in `stats_results*.tsv` (`" ~ Origin$"`), so
renaming it in one place and not the other would break the link between the published
code and the published data. Its values say what it is, and Data S1's column dictionary
defines it.

The three axes of the sampling design are, in both the paper and the code: **fraction**
(VIR, BAC, CUL, LIS, the four prepared samples), **culture dependence** (the two-level
grouping of those four) and **compartment** (the ecological reading, used in the
Introduction and Discussion only).

Some scripts and the knowledge base are in Spanish and are kept as they are, since that
is the version the results came from. `data/mapa_v5_EN.tsv` is a column-for-column
translation for reading; the pipeline reads the Spanish original. Several script headers
are dated change logs and name files that are not here (`Fig5v3.R`, `04_Fig_4_27_8.R`,
"the old pipeline"); those are the record of what changed, not files to look for.

---

## Data flow

Every step reads the previous step's table and writes a new, enriched one.

### vOTU branch

| Script | Reads | Writes |
|---|---|---|
| `votus/01_build_base.py` | representative FASTA, cluster file, CheckV summary | `votu_01_base.tsv` (59,469) |
| `votus/02_add_vibrant.py` | `votu_01_base.tsv`, VIBRANT genome quality | `votu_02_vibrant.tsv` |
| `votus/03_add_coverm.py` | `votu_02_vibrant.tsv`, CoverM TPM table | `votu_03_coverm.tsv` |
| `votus/04_add_phabox_gtdb.py` | `votu_03_coverm.tsv`, PhaBOX summary, GTDB-Tk `bac120` + `ar53` | `votu_04_phabox_gtdb.tsv` |
| `votus/05_add_ct3.py` | `votu_04_phabox_gtdb.tsv`, Cenote-Taker 3 summary | `votu_05_ct3.tsv` |
| `votus/06_filter_catalog.py` | `votu_05_ct3.tsv` | `votu_06_catalog.tsv` (59,469 rows; 11,100 flagged `Passes_catalog`) |

### Gene branch

| Script | Reads | Writes |
|---|---|---|
| `08_defensas_votus.sh` | Prodigal `.faa` | PADLOC CSV + domtblout, DefenseFinder hmmer/genes/systems |
| `genes/01_build_gene_table.py` | Prodigal GFF, VIBRANT annotations | `01_gene_table_L2.tsv` (393,330) |
| `genes/02_add_evidence.py` | `01_gene_table_L2.tsv`, DefenseFinder, PADLOC | `02_gene_table_L3_evidence.tsv` |
| `genes/03_kegg_categories.py` | KEGG pathway/KO tables | `03_ko_kegg_categories.tsv` (253) |
| `genes/04_score_context.py` | `02_gene_table_L3_evidence.tsv` | `03_gene_table_context.tsv` |
| `genes/05_classify.py` | `03_gene_table_context.tsv`, `03_ko_kegg_categories.tsv` | `04_gene_table_classified.tsv` |
| `genes/06_harmonise_categories.py` | `04_gene_table_classified.tsv`, PADLOC `sys/*.yaml` + `hmm_meta.txt` | `06_category_map.tsv` (191 labels) |

### Join, loci, tables and QC

| Script | Reads | Writes |
|---|---|---|
| `07_join_and_report.py` | `votu_06_catalog.tsv`, `04_gene_table_classified.tsv`, `06_category_map.tsv` | `07_long_table_FINAL.tsv`, `07_votu_with_counts.tsv` (59,469) |
| `reconstruir_loci.py` | `07_long_table_FINAL.tsv`, `data/mapa_v5.tsv` | `10_loci_defensa.tsv` (1,718) |
| `08_build_supplementary_tables.py` | `07_votu_with_counts.tsv`, `07_long_table_FINAL.tsv` | `Table_S2_vOTU_catalogue.tsv` (11,100), `Table_S3_AVG_Defence_genes.tsv` (3,405) |
| `98_audit_annotations.py` | `07_long_table_FINAL.tsv` | `98_audit_report.txt` |
| `99_qc.py` | all of the above | `99_qc_report.txt` |

The two tables written by `08_build_supplementary_tables.py` keep the names they had
while the pipeline was written. In the paper they are **Data S1** (vOTU catalogue, 11,100
× 21) and **Data S2** (classified genes, 3,405 × 22), and `10_loci_defensa.tsv` is
**Data S3** (defence loci, 1,718 × 20).

### What feeds R

| R script | Reads | Writes |
|---|---|---|
| `R/01_stats.R` | `07_votu_with_counts.tsv`, `07_long_table_FINAL.tsv` | `stats_results.tsv` |
| `R/02_Fig_2.R` | `07_votu_with_counts.tsv`, `stats_results.tsv` | Figure 2 |
| `R/03_stats_loci.R` | `10_loci_defensa.tsv`, `07_votu_with_counts.tsv`, `06_category_map.tsv` | `07_long_table_loci_for_fig3.tsv`, `votu_n_defence_loci.tsv`, `stats_results_loci_defence.tsv` |
| `R/04_stats_fig3.R` | `07_votu_with_counts.tsv`, `07_long_table_FINAL.tsv`, `10_loci_defensa.tsv` | `stats_results_fig3.tsv`, `fig3_*.tsv`, `stats_report_fig3.txt` |
| `R/05_Fig_3.R` | the four tables above plus the three `stats_results*.tsv` | Figure 3 |
| `R/06_Fig_4.R`, `R/07_Fig_5.R`, `R/08_Fig_6.R` | GenBank and GFF3 files from the folder set in `MAPS_DIR` | Figures 4, 5 and 6 |
| `R/09_supplementary_statistics.R` | Data S1, S2 and S3 as published (`.xlsx` or `.tsv`) | `supplementary_statistics.tsv`, `supplementary_statistics_categories.tsv` |

Figure scripts hard-code no statistic: test results are looked up from the
`stats_results*.tsv` files by section and test name, so rerunning the statistics updates
the figure annotations automatically. `09_supplementary_statistics.R` is independent of
the rest: it recomputes every test in Supplementary Tables S3 and S4 from the three data
files distributed with the paper, and needs `readxl` only for `.xlsx` inputs.

### What feeds the supplementary tables

| Table | Source |
|---|---|
| Table S1 — sequencing and assembly metrics | assembled by hand from the read-processing and QUAST reports |
| Table S2 — category assignment and references | curated; the citations also sit in `config.py` and `data/mapa_v5.tsv` |
| Table S3 — statistics behind Fig. 2 | `R/01_stats.R` and `R/09_supplementary_statistics.R` |
| Table S4 — statistics behind Fig. 3 | `R/09_supplementary_statistics.R` (a, b, c), with the community-level tests from `R/01_stats.R` |
| Table S5 — defence-system families with a complete locus | aggregation of Data S3 |
| Table S6 — evidence behind the gene labels of Fig. 5 | Bakta, VIBRANT, KEGG, Pfam and PHASTEST, compiled by hand |
| Data S1 — vOTU catalogue (11,100) | `08_build_supplementary_tables.py` |
| Data S2 — classified genes (3,405) | `08_build_supplementary_tables.py` |
| Data S3 — defence loci (1,718) | `reconstruir_loci.py` |

`08_build_supplementary_tables.py` checks that the gene table's per-vOTU row counts
reconcile exactly with the catalogue's `N_AMG_genes` and `N_Defence_genes` for all 11,100
vOTUs, and fails loudly if they do not.

---

## 1. vOTU branch — what each virus is

Supplementary Methods S4–S6.

| Script | What it does |
|---|---|
| `01_build_base.py` | Base table from the representative vOTU FASTA (source of truth for IDs). Parses the sampling design from sequence names into four columns (station, fraction group, MDA, mitomycin), traces provirus and MAG provenance from the cluster file, adds CheckV metrics. `Cluster_has_Provirus` is inherited if any cluster member is a CheckV provirus. |
| `02_add_vibrant.py` | VIBRANT lifestyle and quality. Predictions on excised prophages (`<scaffold>_fragment_N`) are reassigned to the parent contig. Where a contig yields several predictions, quality takes the best row and lifestyle takes lysogenic over lytic: a prophage prediction is positive evidence, its absence is not. |
| `03_add_coverm.py` | Abundance (TPM), column names cleaned per station and fraction. Only the two unamplified fractions are used for quantification, but all vOTUs are measured, since libraries are cross-mapped against the whole catalogue. `In_CoverM` separates *absent from the abundance table* from *present but undetected*. |
| `04_add_phabox_gtdb.py` | The four PhaBOX modules (PhaMer, PhaGCN, PhaTYP, CHERRY) plus GTDB-Tk taxonomy for MAG-derived phages, reading both `bac120` and `ar53`. Host assignment follows a hierarchy recorded in `Host_Source`. |
| `05_add_ct3.py` | Cenote-Taker 3 topology, taxonomy and hallmark counts, with merge auditing. Consolidates taxonomy across CT3 and PhaBOX by depth, recording the winner in `Taxonomy_Source`. |
| `06_filter_catalog.py` | The catalogue filter: length with a rescue clause, removal of eukaryotic viruses and virophages, and the count of independent lines of viral evidence, which is recorded but removes no vOTU. Assigns lifestyle in three tiers. |

## 2. Gene branch — what genes it carries

Supplementary Methods S8–S11.

| Script | What it does |
|---|---|
| `01_build_gene_table.py` | ORF-level spine from the Prodigal GFF, with VIBRANT's annotation as separate columns (KO, Pfam, VOG, AMG flag) indexed by position, so each database keeps its own e-value and the AMG flag is read as a field. |
| `02_add_evidence.py` | DefenseFinder and PADLOC results as evidence columns. Nothing is classified here: each tool keeps its own opinion in its own columns. PADLOC hits use PADLOC's per-profile thresholds from `hmm_meta.txt`. Also flags whether each hit belongs to a complete multi-gene system. |
| `03_kegg_categories.py` | KO → KEGG category lookup for AMGs. A KO may belong to several categories; the primary one is where it occurs in the greatest number of distinct pathways, with a declared tie-break. |
| `04_score_context.py` | Genomic-context auxiliary score (viral hallmark genes within ±10 ORFs), computed for every ORF so AMGs and defence genes meet the same criterion. |
| `05_classify.py` | The single classification decision. Produces `gene_class`, `class_source` and the raw category label. |
| `06_harmonise_categories.py` | Resolves each raw label to one of six mechanistic categories, trying the PADLOC system name, then the profile's protein names, then the label. Writes `06_category_map.tsv` with the route taken for every label. |

### The classification rule

One precedence, all of it in `05_classify.py`:

| Rule | Source | Result |
|---|---|---|
| 0 | core-function veto | evaluated once per protein, before any rule |
| 1 | DefenseFinder or PADLOC assigned a defence system | `Defence` |
| 2 | KO name or Pfam name matches `SEARCH_TERMS` | `Defence` |
| 3 | VIBRANT AMG flag + KEGG lookup | `AVG` (reported as AMG) |
| 4 | — | `Other` |

Patterns are matched only against KO and Pfam names, never against VOG names, which
describe a reference protein from another organism. Rule 2 exists because the specialist
tools model bacterial anti-phage systems and do not search for lysogenic conversion,
superinfection exclusion, phage communication or SOS modulation. Rule 1 outranks rule 3
because a profile built on a described system is stronger evidence than a generic
auxiliary-gene flag; in this dataset that reclassifies 579 genes VIBRANT had flagged as
metabolic, 94% of them DNA methyltransferases a specialist places inside a
restriction–modification system.

**Core-function veto.** A defence assignment is rejected when the protein's own
annotation identifies it as a core replication, repair or metabolic enzyme, or a prophage
lifecycle regulator — DNA polymerase A, primases, SNF2 and RecB-family helicases and
preQ₀ enzymes all appear under defence labels. It is evaluated once, on the protein,
before any rule dispatches. A vetoed protein falls through to rule 3 and is counted as an
AMG. All 319 vetoes are recorded in `class_source` and exported in Data S2.

### Defence categories

Six mechanism-based categories following Georjon & Bernheim (2023) *Nat Rev Microbiol*
21:686–700: `CRISPR-Cas` · `Nucleic acid degradation` · `Genome modification` ·
`Abortive infection` · `Other / uncharacterised systems` ·
`Host manipulation (phage-encoded)`.

A system is assigned a mechanism only where a primary source describes its molecular
activity; the rest go to *Other / uncharacterised systems*, which therefore reflects the
state of the literature on the family and not the quality of the gene's annotation.
Assignments were checked against Doron et al. (2018), Hampton, Watson & Fineran (2020)
and Georjon & Bernheim (2023), plus later sources where a mechanism was resolved
afterwards. The citation for each sits beside it in `config.py` and `data/mapa_v5.tsv`.

Categories are **system-level, not gene-level**: the methyltransferase subunit of a
restriction–modification system appears under *Nucleic acid degradation*. The exception
is a call the locus reconstruction rescinds.

BREX, DISARM and MADS are placed under *Genome modification*, following the established
methylation-based part of their mechanism. MADS is counted at gene level only and
excluded from locus reconstruction, because no vOTU carries enough of its components for
the architecture to be evaluated (37 genes in 36 vOTUs).

> **One correction.** The MADS rows in `data/mapa_v5.tsv` describe the family as a
> probable annotation artefact and cite a paper on plant MADS-box transcription factors.
> MADS is in fact a described eight-gene anti-phage system (*Cell Host Microbe* 2024,
> doi:10.1016/j.chom.2024.07.005). The decision those rows record is the one the paper
> applies and is unchanged; only the stated reason was wrong.

---

## 3. Defence loci — `reconstruir_loci.py`

Supplementary Methods S12.

Groups defence genes into loci by contig position and tests each against the architecture
its family requires. The architectures are the machine-readable expressions in the `ARQ`
dictionary at the top of the script; the curated evidence per family — window, level of
evidence, confidence, what to do with a fragment — is `data/mapa_v5.tsv`, 22 columns over
203 families.

A locus counts only when every required component is present: GajA and GajB for Gabija, a
reverse transcriptase in the window for a retron, a methyltransferase and an endonuclease
for a type II RM system. **Each required component must come from a distinct gene** —
otherwise a label like `Toxin_Antitoxin_and_Abi` satisfies both the toxin and the
antitoxin term by itself. Positive terms are matched to distinct genes by bipartite
matching; negative terms are evaluated over the whole block, since an absence has to hold
across the neighbourhood.

**Where the evidence rescinds a call, the category follows.** A methyltransferase with no
endonuclease is not an RM system: it becomes an orphan methyltransferase under *Genome
modification*.

Each locus also carries an **axis**, named from the phage's perspective:
*self-protection*, for elements compatible with protecting the phage genome, or
*exclusion*, for systems acting against incoming genetic material. Loci that are neither
— lifecycle regulation, host manipulation, unresolved candidates — keep their own axis
and are reported in Data S3 but excluded from Figure 3. A complete locus is **confirmed**
when its system is experimentally characterised with a diagnostic architecture, and a
**candidate** otherwise.

---

## 4. Join, tables and QC

| Script | What it does |
|---|---|
| `07_join_and_report.py` | Joins both branches and applies the category map. Per-vOTU gene counts are derived here, once, by aggregation. Filtering is left to reporting time: the long table carries `Passes_catalog`, `Conservative_10kb` and `passes_context` as flags, so every sensitivity analysis is a `filter()` in R. |
| `08_build_supplementary_tables.py` | The vOTU table (Data S1) and the gene table (Data S2), checking that the two reconcile. Applies the AMG and Virulent presentation labels. |
| `98_audit_annotations.py` | Audits every gene against its own annotation: no defence gene carrying a core-function identity, flags any whose annotation contains replication, repair, ribosomal or structural vocabulary. |
| `99_qc.py` | Row counts preserved, no duplicate IDs, derived counts equal to their aggregation, flags correctly nested, category map complete — plus the manuscript numbers in one block. Exits non-zero on failure. |

---

## 5. Statistics and figures — `R/`

Supplementary Methods S14–S15.

| Script | What it does |
|---|---|
| `01_stats.R` | Catalogue description, lifestyle and its controls, alpha and beta diversity, biogeography, gene content and abundance. PERMANOVA p-values are exact: with four libraries per group all 35 partitions are enumerated, so the smallest attainable p is 1/35 = 0.029. |
| `02_Fig_2.R` | Figure 2 — vOTUs and abundance by fraction, Bray–Curtis PCoA and PERMANOVA, lifestyle by culture dependence, station and fraction overlap, taxonomy and predicted host. |
| `03_stats_loci.R` | Defence statistics at **locus** level and the category overrides that follow the reconstruction. Skipping it silently reverts Figure 3 to counting individual genes. |
| `04_stats_fig3.R` | Carrier prevalence, the AMG side of the category test, axis composition across fractions and the demography of the two axes. |
| `05_Fig_3.R` | Figure 3, five panels. Looks up every statistic in the `stats_results*.tsv` files and refuses to draw crude odds ratios. |
| `06_Fig_4.R`, `07_Fig_5.R`, `08_Fig_6.R` | The genomic maps. Gene categories are assigned by hand from the classification, the locus and the annotation, and are documented gene by gene inside each script. |
| `09_supplementary_statistics.R` | Recomputes Supplementary Tables S3 and S4 from Data S1–S3 alone, in base R + MASS. |
| `palette_paper.R` | Colours, category names and factor orders for Figures 2 and 3. |
| `palette_maps.R` | The same for Figures 4, 5 and 6. One file for the three, so a colour means the same thing in all of them; each script prints which palette it loaded. |

Odds ratios for the category-level contrasts are **Cochran–Mantel–Haenszel**: culture
dependence and lifestyle are confounded in this design, so each contrast is estimated
with the other held constant. The crude Fisher values are written alongside, in their own
sections of the results files, and are not what the figures draw.

Genomic-map figures use representative genomes from the catalogue. Their boundaries,
lengths, quality and taxonomy come from the pipeline, but ORF coordinates and gene labels
are annotated separately with Cenote-Taker 3 or Bakta and cross-checked against VIBRANT
and PHASTEST, so labels there are finer than the six mechanistic categories.

Each figure script writes PDF, EPS, TIFF (600 dpi, LZW) and PNG into an `ISME_figs/`
folder under the working directory, named `fig2` to `fig6` to match how they are cited.
The PDFs are fully vector with the fonts embedded; the EPS of Figure 2 rasterises the two
semi-transparent boxplot panels, because EPS cannot carry transparency.

---

## Running it

This is the order the scripts were run in. Without the raw inputs it will not run from
here.

```bash
export PYTHONPATH=/path/to/pipeline
python3 votus/01_build_base.py          # … through votus/06
python3 genes/01_build_gene_table.py    # … through genes/06
python3 07_join_and_report.py
python3 08_build_supplementary_tables.py
python3 reconstruir_loci.py
python3 98_audit_annotations.py
python3 99_qc.py
```

`08_defensas_votus.sh` runs separately, before the gene branch: it produces the PADLOC and
DefenseFinder calls that `genes/02_add_evidence.py` reads.

Then, with the working directory set to the folder holding the pipeline tables:

```r
source("01_stats.R")        # -> stats_results.tsv
source("02_Fig_2.R")        # -> Figure 2
source("03_stats_loci.R")   # -> stats_results_loci_defence.tsv
source("04_stats_fig3.R")   # -> stats_results_fig3.tsv
source("05_Fig_3.R")        # -> Figure 3
source("06_Fig_4.R"); source("07_Fig_5.R"); source("08_Fig_6.R")
source("09_supplementary_statistics.R")
```

The R scripts take their input names with or without the `_28_8` suffix the tables were
exported with, and print which file they opened. The three map scripts read GenBank and
GFF3 files from `MAPS_DIR`, set at the top of each, and the palette from `PALETTE`; those
two constants are the only paths to set by hand.

Every script prints input and output counts, filtering statistics and validation checks,
so each stage can be verified on its own. In the gene branch the script numbers and the
output-file numbers do not line up (`genes/04` writes `03_gene_table_context.tsv`); the
authoritative mapping is the `F_*` constants in `config.py`.

---

## Inputs

Not included, due to size. Expected in the locations set at the top of `config.py`, all of
which can be overridden with the `WAP_BASE`, `WAP_CCBB` and `PADLOC_DATA` environment
variables. Every header is asserted on read, so a shifted column fails loudly instead of
producing plausible wrong numbers.

| Directory | Contents |
|---|---|
| `input/` | representative vOTU FASTA, cluster file |
| `checkv/` | CheckV quality summary |
| `vibrant/` | VIBRANT genome quality, annotations, AMG pathways, prophage coordinates |
| `coverm/` | CoverM TPM abundance table |
| `phabox/` | PhaBOX final prediction summary |
| `ct3/` | Cenote-Taker 3 virus summary |
| `gtdbtk/` | GTDB-Tk `bac120` and `ar53` summaries |
| `defence/` | DefenseFinder HMMER/systems/genes tables, PADLOC CSV and domtblout |

The GenBank and GFF3 files behind the genomic maps are not included either; they are read
from whatever folder `MAPS_DIR` points at.

## Outputs

| File | Contents | In the paper |
|---|---|---|
| `votu_06_catalog.tsv` | all 59,469 representatives, flagged; 11,100 pass | — |
| `07_long_table_FINAL.tsv` | one row per AMG/defence gene, with vOTU metadata and final category | — |
| `07_votu_with_counts.tsv` | the vOTU table plus `n_orfs`, `n_AVG`, `n_Defence` | — |
| `10_loci_defensa.tsv` | 1,718 defence loci, with axis, completeness and any reassignment | Data S3 |
| `Table_S2_vOTU_catalogue.tsv` | 11,100 vOTUs | Data S1 |
| `Table_S3_AVG_Defence_genes.tsv` | 3,405 genes (1,441 AMG + 1,964 defence) | Data S2 |
| `06_category_map.tsv` | raw defence label → mechanistic category, with the route | — |
| `stats_results*.tsv` | every test behind Figures 2 and 3 | Tables S3, S4 |
| `98_audit_report.txt` | annotation audit | — |
| `99_qc_report.txt` | consistency checks and the manuscript numbers | — |

---

## Software and parameters

CheckV · VIBRANT · PhaBOX · CoverM · GTDB-Tk · Cenote-Taker 3 · Prodigal · DefenseFinder ·
PADLOC · Bakta · PHASTEST · vConTACT3 · Python 3 (pandas, scipy) · R (tidyverse, vegan,
MASS, cowplot, scales, ggupset, gggenes, readxl). Versions in `requirements.txt` and
`environment.txt`.

| Parameter | Value | Source |
|---|---|---|
| HMM e-value cutoff | 1e-5 | as VIBRANT (Kieft, Zhou & Anantharaman 2020, *Microbiome* 8:90) |
| genomic context window | ±10 ORFs | criterion from DRAM-v (Shaffer et al. 2020, *NAR* 48:8883) |
| auxiliary score retained | ≤ 3 | `DRAM-v distill` default |
| minimum contig length | 5,000 bp | flag, never a filter |
| conservative subset | ≥ 10,000 bp or complete | Pratama et al. 2021, *PeerJ* 9:e11447 |
| vOTU clustering | ≥ 95% ANI over ≥ 85% of the shorter sequence | standard practice |

DRAM-v was not run: its auxiliary score is a positional rule, so the criterion and its
default threshold were computed directly from the ORF table in `genes/04_score_context.py`.
The score is recorded but not applied as a filter, and is exported per gene so readers can
apply it if they wish.

---

## Design decisions worth knowing

- **Nothing is deleted.** No script removes rows; filters are flags, applied at reporting
  time.
- **Physical evidence outranks inference.** A co-assembled MAG over a database host
  prediction, a CheckV provirus boundary over a lifestyle prediction, a specialist call
  over a name match, a complete locus over an isolated component. The source of every call
  is recorded.
- **Assertions abort.** Fragment ORF numbering is verified to fill its declared parent
  range exactly, because shifted identifiers would pass a mere existence check while
  pointing at the wrong gene.
- Per-vOTU gene counts are derived once and never recomputed from a different source.
- The number of independent lines of viral evidence is descriptive: reported per vOTU,
  removes none.
- p-values below 1e-15 are reported from the test statistic rather than from R's
  floating-point floor of 2.22e-16.

---

## Licence and citation

Code: MIT. Curated data files (`data/mapa_v5.tsv`, `data/mapa_v5_EN.tsv`): CC BY 4.0. See
`LICENSE`.

This repository is archived on Zenodo; cite the DOI there rather than the GitHub URL,
since the DOI points at a fixed version. Data S1–S3 are distributed with the paper. The
raw reads are deposited at the European Nucleotide Archive under study accession
PRJEB128412 and are released on publication.
