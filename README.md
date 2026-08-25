# WAP-virosphere

Pipeline for the construction, curation, functional annotation and
statistical analysis of a viral operational taxonomic unit (vOTU)
catalogue from the Western Antarctic Peninsula, and the code behind every
number and figure in the associated paper.

Raw inputs are not included, so the pipeline is not runnable end to end
from this repository alone.

## The paper

Polar viromics has focused almost exclusively on free-living viral
particles, leaving the integrated prophage reservoir largely unexplored.
We combined four fractions — free-living viral concentrates, environmental
bacterial fractions, cultivated bacterial communities, and mitomycin
C-induced lysates — with hybrid short- and long-read assembly across four
stations, yielding 11,100 non-redundant vOTUs. Viral assemblages
partitioned sharply by fraction, and so did their functional content.
Synteny mapping of co-resident prophages in *Vibrio diabolicus* and
*V. alginolyticus* MAGs revealed stacked restriction–modification,
addiction and diversity-generating retroelement systems.

### Key numbers

From `99_qc.py`. The manuscript copies these rather than re-deriving them.

| | |
|---|---|
| vOTUs assembled | 59,469 |
| vOTUs in the final catalogue | 11,100 |
| supported by ≥ 2 lines of viral evidence | 10,671 (96.1%) |
| conservative subset (≥ 10 kb or circular) | 3,810 |
| auxiliary metabolic genes (AMG) | 1,441 |
| defence-associated genes | 1,964 |
| reconstructed defence loci | 1,718 |
| complete defence loci | 808 |
| complete loci on the two defence axes | 773 |

Defence genes by category, gene level:

| Category | n | % |
|---|---|---|
| Nucleic acid degradation | 840 | 42.8 |
| Abortive infection | 383 | 19.5 |
| Other / uncharacterised systems | 386 | 19.7 |
| Genome modification | 318 | 16.2 |
| Host manipulation (phage-encoded) | 31 | 1.6 |
| CRISPR-Cas | 6 | 0.3 |

Complete defence loci by axis:

| Axis | n | % |
|---|---|---|
| Counter-defence | 599 | 77.5 |
| Defence | 174 | 22.5 |

Gene-level and locus-level distributions differ because a locus whose
restriction–modification call is rescinded for lack of an endonuclease is
recounted as an orphan methyltransferase under *Genome modification*.

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
08_build_supplementary_tables.py   Tables S2 and S3
08_defensas_votus.sh               PADLOC + DefenseFinder screening
98_audit_annotations.py            every gene against its own annotation
99_qc.py                           consistency checks + manuscript numbers
data/
  mapa_v5.tsv                      curated family knowledge base
  mapa_v5_EN.tsv                   English translation, for reading only
R/                                 statistics and figures
requirements.txt  environment.txt  README.md  LICENSE
```

All paths, thresholds and category definitions live in `config.py`. If any
of them appears elsewhere, that is a bug.

**Nomenclature.** The manuscript and the supplementary tables say **AMG**
(auxiliary metabolic gene). Inside the code the same set of genes is called
`AVG` (auxiliary viral gene), which is what the variable was named while
the pipeline was written: `gene_class == "AVG"`, `n_AVG`, and the output
filename `Table_S3_AVG_Defence_genes.tsv`. The relabel to AMG is applied
once, in `08_build_supplementary_tables.py`, at the boundary between the
pipeline and the reader. The internal names were left alone so that the
code here is the code the results came from.

The same applies to lifestyle. The tables report **Virulent**, following
the convention that a *phage* is virulent or temperate while a *cycle* is
lytic or lysogenic; the internal value is `Lytic`, which is what VIBRANT
and PhaBOX return.

Some scripts and the knowledge base were written in Spanish and are kept
as they are, since that is the version the results came from.
`data/mapa_v5_EN.tsv` is a column-for-column translation for reading; the
pipeline reads the Spanish original.

**On the comments.** Many comments record why a rule exists, including the
counts from the audit that produced it and the cost of the decision. They
are part of the documentation, not leftovers.

---

## Data flow

Every step reads the previous step's table and writes a new, enriched one.
Row counts are given where they are informative.

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
| `07_join_and_report.py` | `votu_06_catalog.tsv`, `04_gene_table_classified.tsv`, `06_category_map.tsv` | `07_long_table_FINAL.tsv` (7,364), `07_votu_with_counts.tsv` (59,469) |
| `reconstruir_loci.py` | `07_long_table_FINAL.tsv`, `data/mapa_v5.tsv` | `10_loci_defensa.tsv` (1,718) |
| `08_build_supplementary_tables.py` | `07_votu_with_counts.tsv`, `07_long_table_FINAL.tsv` | `Table_S2_vOTU_catalogue.tsv` (11,100), `Table_S3_AVG_Defence_genes.tsv` (3,405) |
| `98_audit_annotations.py` | `07_long_table_FINAL.tsv` | `98_audit_report.txt` |
| `99_qc.py` | all of the above | `99_qc_report.txt` |

### What feeds R

| R script | Reads | Writes |
|---|---|---|
| `R/01_stats_F.R` | `07_votu_with_counts.tsv`, `07_long_table_FINAL.tsv` | `stats_results.tsv` |
| `R/02_stats_loci_defence_F.R` | `10_loci_defensa.tsv`, `07_votu_with_counts.tsv`, `06_category_map.tsv` | `07_long_table_loci_for_fig3.tsv`, `votu_n_defence_loci.tsv`, `stats_results_loci_defence.tsv` |
| `R/03_Fig2.R` | `07_votu_with_counts.tsv`, `stats_results.tsv` | Figure 2 (PNG, TIFF, EPS) |

Figure scripts hard-code no statistic: test results are looked up from
`stats_results.tsv` by section and test name, so rerunning the statistics
updates the figure annotations automatically.

### What feeds the supplementary tables

| Table | Source |
|---|---|
| S1 — sequencing and assembly metrics | assembled by hand from the read-processing and QUAST reports |
| S2 — vOTU catalogue | `08_build_supplementary_tables.py`, from `07_votu_with_counts.tsv` |
| S3 — AMG and defence genes | `08_build_supplementary_tables.py`, from `07_long_table_FINAL.tsv` |
| S4 — category assignment and references | curated; the citations also sit in `config.py` and `data/mapa_v5.tsv` |
| S5 — defence loci | `10_loci_defensa.tsv` |
| S6 — loci summarised per family | aggregation of S5 |
| S7 — statistical results | `stats_results.tsv` from `R/01_stats_F.R` |
| S8 — gene annotations for the genomic maps | Cenote-Taker 3 / Bakta, cross-checked against VIBRANT and PHASTEST |

`08_build_supplementary_tables.py` checks that S3's per-vOTU row counts
reconcile exactly with S2's `N_AMG_genes` and `N_Defence_genes` for all
11,100 vOTUs, and fails loudly if they do not.

---

## 1. vOTU branch — what each virus is

| Script | What it does |
|---|---|
| `01_build_base.py` | Base table from the representative vOTU FASTA (source of truth for IDs). Parses the sampling design from sequence names into four columns (zone, fraction group, MDA, mitomycin), traces provirus and MAG origin from the cluster file, adds CheckV metrics. `Cluster_has_Provirus` is inherited if any cluster member is a CheckV provirus. |
| `02_add_vibrant.py` | VIBRANT lifestyle, quality and `Vibrant_predicted`. Predictions on excised prophages (`<scaffold>_fragment_N`) are reassigned to the parent contig. Where a contig yields several predictions, quality takes the best row and lifestyle takes lysogenic over lytic: a prophage prediction is positive evidence, its absence is not. |
| `03_add_coverm.py` | Abundance (TPM), column names cleaned per zone and fraction. Only the two non-MDA libraries are used for quantification, but all vOTUs are measured, since libraries are cross-mapped against the whole catalogue. `In_CoverM` separates *absent from the abundance table* from *present but undetected*. |
| `04_add_phabox_gtdb.py` | The four PhaBOX modules (PhaMer: viral; PhaGCN: taxonomy; PhaTYP: lifestyle; CHERRY: host) plus GTDB-Tk taxonomy for MAG-derived phages, reading both `bac120` and `ar53`. Host assignment follows a hierarchy recorded in `Host_Source`. |
| `05_add_ct3.py` | Cenote-Taker 3 topology, taxonomy and hallmark counts, with merge auditing. Consolidates taxonomy across CT3 and PhaBOX by depth, recording the winner in `Taxonomy_Source`. |
| `06_filter_catalog.py` | The catalogue filter: viral evidence, length with a rescue clause, removal of eukaryotic viruses and virophages. Assigns lifestyle in three tiers. Prints a cascading filtering summary. |

## 2. Gene branch — what genes it carries

| Script | What it does |
|---|---|
| `01_build_gene_table.py` | ORF-level spine from the Prodigal GFF, with VIBRANT's annotation as separate columns (KO, Pfam, VOG, AMG flag) indexed by position, so each database keeps its own e-value and the AMG flag is read as a field. |
| `02_add_evidence.py` | DefenseFinder and PADLOC results as evidence columns. Nothing is classified here: each tool keeps its own opinion in its own columns. PADLOC hits use PADLOC's per-profile thresholds from `hmm_meta.txt`. Also flags whether each hit belongs to a complete multi-gene system. |
| `03_kegg_categories.py` | KO → KEGG category lookup for AMGs. A KO may belong to several categories; the primary one is where it occurs in the greatest number of distinct pathways, with a declared tie-break. Full membership is kept alongside. |
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

Patterns are matched only against KO and Pfam names, never against VOG
names, which describe a reference protein from another organism. Rule 2
exists because the specialist tools model bacterial anti-phage systems and
do not search for lysogenic conversion, superinfection exclusion, phage
communication or SOS modulation.

Rule 1 outranks rule 3 because a profile built on a described system is
stronger evidence than a generic auxiliary-gene flag. In this dataset that
reclassifies 579 genes that VIBRANT had flagged as metabolic, 94% of them
DNA methyltransferases that a specialist places inside a
restriction–modification system.

**Core-function veto.** A defence assignment is rejected when the
protein's own annotation identifies it as a core replication, repair or
metabolic enzyme, or a prophage lifecycle regulator — DNA polymerase A,
primases, SNF2 and RecB-family helicases and preQ₀ enzymes all appear
under defence labels. It is evaluated once, on the protein, before any
rule dispatches, so that a protein with no specialist call cannot reach
rule 2 unchecked. A vetoed protein falls through to rule 3 and is counted
as an AMG. Every veto is recorded in `class_source` and exported in
Table S3.

### Defence categories

Six mechanism-based categories following Georjon & Bernheim (2023)
*Nat Rev Microbiol* 21:686–700: `CRISPR-Cas` · `Nucleic acid degradation` ·
`Genome modification` · `Abortive infection` ·
`Other / uncharacterised systems` · `Host manipulation (phage-encoded)`.

A system is assigned a mechanism only where a primary source describes its
molecular activity; the rest go to *Other*. Assignments were checked
against Doron et al. (2018), Hampton, Watson & Fineran (2020) and Georjon
& Bernheim (2023), plus later sources where a mechanism was resolved
afterwards. The citation for each sits beside it in `config.py` and
`data/mapa_v5.tsv`, and is collected in Table S4.

Categories are **system-level, not gene-level**: the methyltransferase
subunit of a restriction–modification system appears under *Nucleic acid
degradation*. The exception is a call the locus reconstruction rescinds.

---

## 3. Defence loci — `reconstruir_loci.py`

Groups defence genes into loci by contig position and tests each against
the architecture its family requires. The architectures are the
machine-readable expressions in the `ARQ` dictionary at the top of the
script; the curated evidence per family — window, level of evidence,
confidence, what to do with a fragment — is `data/mapa_v5.tsv`, 22 columns
over 203 families.

A locus counts only when every required component is present: GajA and
GajB for Gabija, a reverse transcriptase in the window for a retron, a
methyltransferase and an endonuclease for a type II RM system. **Each
required component must come from a distinct gene** — otherwise a label
like `Toxin_Antitoxin_and_Abi` satisfies both the toxin and the antitoxin
term by itself. Positive terms are matched to distinct genes by bipartite
matching; negative terms are evaluated over the whole block, since an
absence has to hold across the neighbourhood.

**Where the evidence rescinds a call, the category follows.** A
methyltransferase with no endonuclease is not an RM system: it becomes an
orphan methyltransferase under *Genome modification*.

Each locus also carries an **axis**: counter-defence, protecting the phage
genome from host systems, or defence proper, acting against incoming
genetic material.

---

## 4. Join, tables and QC

| Script | What it does |
|---|---|
| `07_join_and_report.py` | Joins both branches and applies the category map. Per-vOTU gene counts are derived here, once, by aggregation. Filtering is left to reporting time: the long table carries `Passes_catalog`, `Conservative_10kb` and `passes_context` as flags, so every sensitivity analysis is a `filter()` in R. |
| `08_build_supplementary_tables.py` | Tables S2 (one row per vOTU) and S3 (one row per AMG/defence gene), checking that the two reconcile. Applies the AMG and Virulent presentation labels. |
| `98_audit_annotations.py` | Audits every gene against its own annotation: no defence gene carrying a core-function identity, flags any whose annotation contains replication, repair, ribosomal or structural vocabulary, prints the full inventory of both sets. |
| `99_qc.py` | Row counts preserved, no duplicate IDs, derived counts equal to their aggregation, flags correctly nested, category map complete — plus the manuscript numbers in one block. Exits non-zero on failure. |

---

## 5. Figures — `R/`

| Script | What it does |
|---|---|
| `01_stats_F.R` | Alpha and beta diversity, Fisher, Wilcoxon and χ² with Benjamini–Hochberg, biogeography, completeness-versus-lifestyle bias control. Writes `stats_results.tsv`. |
| `02_stats_loci_defence_F.R` | Defence statistics at **locus** level and the category reassignment. Runs between `reconstruir_loci.py` and the Figure 3 script; skipping it silently reverts Figure 3 to counting individual genes. |
| `03_Fig2.R` | Fig. 2 — richness and abundance by fraction, Bray–Curtis PCoA and PERMANOVA, lifestyle, zone and fraction overlap, taxonomy, predicted host. |

Genomic-map figures use representative genomes from the catalogue. Their
boundaries, lengths, quality and taxonomy come from the pipeline, but ORF
coordinates and gene labels are annotated separately with Cenote-Taker 3
or Bakta and cross-checked against VIBRANT and PHASTEST, so labels there
are finer than the six mechanistic categories.

All figure scripts export PNG, TIFF (600 dpi, LZW) and EPS.

---

## Running it

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

`08_defensas_votus.sh` runs separately, before the gene branch: it
produces the PADLOC and DefenseFinder calls that `genes/02_add_evidence.py`
reads.

Every script prints input and output counts, filtering statistics and
validation checks, so each stage can be verified on its own.

In the gene branch the script numbers and the output-file numbers do not
line up (`genes/04` writes `03_gene_table_context.tsv`, `genes/05` writes
`04_gene_table_classified.tsv`). The authoritative mapping is the `F_*`
constants in `config.py`.

---

## Inputs

Not included, due to size. Expected in the locations set at the top of
`config.py`. Every header is asserted on read, so a shifted column fails
loudly instead of producing plausible wrong numbers.

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
| `data/MAPAS/` | GenBank and GFF3 for the genomic map figures |

## Outputs

| File | Contents |
|---|---|
| `votu_06_catalog.tsv` | all 59,469 representatives, flagged; 11,100 pass |
| `07_long_table_FINAL.tsv` | one row per AMG/defence gene, with vOTU metadata and final category |
| `07_votu_with_counts.tsv` | the vOTU table plus `n_orfs`, `n_AVG`, `n_Defence` |
| `10_loci_defensa.tsv` | 1,718 defence loci, with axis, completeness and any reassignment |
| `Table_S2_vOTU_catalogue.tsv` | Table S2 — 11,100 vOTUs, with `N_AMG_genes` and `N_Defence_genes` |
| `Table_S3_AVG_Defence_genes.tsv` | Table S3 — 3,405 genes (1,441 AMG + 1,964 defence) |
| `06_category_map.tsv` | raw defence label → mechanistic category, with the route |
| `98_audit_report.txt` | annotation audit |
| `99_qc_report.txt` | consistency checks and the manuscript numbers |

---

## Software and parameters

CheckV · VIBRANT · PhaBOX · CoverM · GTDB-Tk · Cenote-Taker 3 · Prodigal ·
DefenseFinder · PADLOC · Bakta · PHASTEST · vConTACT3 · Python 3 (pandas,
scipy) · R (tidyverse, vegan, MASS, gggenes, pheatmap, ggupset).
Versions in `requirements.txt` and `environment.txt`.

| Parameter | Value | Source |
|---|---|---|
| HMM e-value cutoff | 1e-5 | as VIBRANT (Kieft, Zhou & Anantharaman 2020, *Microbiome* 8:90) |
| genomic context window | ±10 ORFs | criterion from DRAM-v (Shaffer et al. 2020, *NAR* 48:8883) |
| auxiliary score retained | ≤ 3 | `DRAM-v distill` default |
| minimum contig length | 5,000 bp | flag, never a filter |
| conservative subset | ≥ 10,000 bp or circular | Pratama et al. 2021, *PeerJ* 9:e11447 |
| vOTU clustering | ≥ 95% ANI over ≥ 85% of the shorter sequence | standard practice |

DRAM-v was not run: its auxiliary score is a positional rule, so the
criterion and its default threshold were computed directly from the ORF
table in `genes/04_score_context.py`. The score is recorded but not
applied as a filter — its pass rate depends on contig length, which
differs between the fractions being compared — and is exported per gene so
readers can apply it if they wish.

---

## Design decisions worth knowing

- **Nothing is deleted.** No script removes rows; filters are flags,
  applied at reporting time.
- **Physical evidence outranks inference.** A co-assembled MAG over a
  database host prediction, a CheckV provirus boundary over a lifestyle
  prediction, a specialist call over a name match, a complete locus over
  an isolated component. The source of every call is recorded.
- **Assertions abort.** Fragment ORF numbering is verified to fill its
  declared parent range exactly, because shifted identifiers would pass a
  mere existence check while pointing at the wrong gene.
- Per-vOTU gene counts are derived once and never recomputed from a
  different source.
- Paths in the R scripts are absolute and reflect the machine they ran on.

## Data availability

Raw reads: [BioProject accession]
Assemblies, catalogue and annotation tables: [repository, DOI]

## Citation

Please cite the associated paper (details on publication). Code archived
at Zenodo: DOI [xxx].

## License

[MIT / GPL-3.0 / CC-BY-4.0 — choose one and match the LICENSE file]
