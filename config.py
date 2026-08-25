#!/usr/bin/env python3
"""
config.py — single source of truth for BOTH branches of the pipeline.

    votus/   one row per vOTU   (catalogue, sampling design, abundance)
    genes/   one row per ORF    (annotation, evidence, context, class)

If a path, a threshold or a category definition appears anywhere else in
the pipeline, that is a bug.

Usage:
    export PYTHONPATH=~/gabi/vOTUs_total/pipeline_v2
    python3 votus/01_build_base.py
"""

import os

# ==========================================================================
# 1. BASE DIRECTORIES
# ==========================================================================
# WAP_BASE  : directorio con los datos de entrada (ver README)
# PADLOC_DATA: padloc-db v2.0.0 (trae hmm_meta.txt y sys/)
BASE = os.environ.get("WAP_BASE", os.path.expanduser("~/gabi/vOTUs_total"))
CCBB = os.environ.get("WAP_CCBB", os.path.expanduser("~/gabi/ccbb"))
PADLOC_DATA = os.environ.get(
    "PADLOC_DATA", os.path.expanduser("~/miniconda3/envs/padloc/data"))

VIBRANT_RES = os.path.join(
    BASE, "vibrant_results/VIBRANT_vOTUs_catalogue_rep_seq",
    "VIBRANT_results_vOTUs_catalogue_rep_seq")
DEFENCE_MERGED = os.path.join(BASE, "analisis_defensas_votus/merged")
OUT_DIR = os.path.join(BASE, "pipeline_v2/out")

os.makedirs(OUT_DIR, exist_ok=True)

# ==========================================================================
# 2. INPUTS — vOTU branch
# ==========================================================================
F_CATALOGUE_FNA = os.path.join(BASE, "vOTUs_catalogue_rep_seq.fasta")
F_CLUSTER_TSV = os.path.join(BASE, "vOTUs_catalogue_cluster.tsv")
F_CHECKV = os.path.join(BASE, "master_quality_summary.tsv")
F_COVERM = os.path.join(BASE, "abundance_coverm/master_vOTUs_tpm_table.tsv")
F_PHABOX = os.path.join(BASE, "phabox/final_prediction_summary.tsv")
F_CT3 = os.path.join(
    BASE, "ct3_results/vOTU_analysis_results/vOTU_viral",
    "vOTU_viral_virus_summary.tsv")
F_VIB_QUALITY = os.path.join(
    VIBRANT_RES, "VIBRANT_genome_quality_vOTUs_catalogue_rep_seq.tsv")

# GTDB-Tk ships bacteria and archaea in SEPARATE files. The old pipeline
# read a single 'gtdbtk.summary.tsv', so any archaeal MAG silently lost its
# taxonomy. Both are read and concatenated in votus/04.
F_GTDB_BAC = os.path.join(CCBB, "gtdbtk_results/gtdbtk.bac120.summary.tsv")
F_GTDB_AR = os.path.join(CCBB, "gtdbtk_results/gtdbtk.ar53.summary.tsv")

# ==========================================================================
# 3. INPUTS — gene branch
# ==========================================================================
F_PRODIGAL_GFF = os.path.join(DEFENCE_MERGED, "vOTUs_prodigal.gff")
F_PRODIGAL_FAA = os.path.join(DEFENCE_MERGED, "vOTUs_prodigal.faa")

F_VIB_ANNOT = os.path.join(
    VIBRANT_RES, "VIBRANT_annotations_vOTUs_catalogue_rep_seq.tsv")
F_VIB_AMG_PATH = os.path.join(
    VIBRANT_RES, "VIBRANT_AMG_pathways_vOTUs_catalogue_rep_seq.tsv")
F_VIB_PROPHAGE = os.path.join(
    VIBRANT_RES,
    "VIBRANT_integrated_prophage_coordinates_vOTUs_catalogue_rep_seq.tsv")

F_DF_HMMER = os.path.join(DEFENCE_MERGED, "defense_finder_hmmer_all.tsv")
F_DF_SYSTEMS = os.path.join(DEFENCE_MERGED, "defense_finder_systems_all.tsv")
F_DF_GENES = os.path.join(DEFENCE_MERGED, "defense_finder_genes_all.tsv")
F_PADLOC_CSV = os.path.join(DEFENCE_MERGED, "padloc_all.csv")
F_PADLOC_DOMTBL = os.path.join(DEFENCE_MERGED, "padloc_domtblout_all.txt")
F_PADLOC_HMM_META = os.path.join(PADLOC_DATA, "hmm_meta.txt")
F_PADLOC_SYS_DIR = os.path.join(PADLOC_DATA, "sys")

# ==========================================================================
# 3b. EXPECTED SCHEMAS  (asserted on read — fail loud, never shift silently)
# ==========================================================================
# VIBRANT annotations: 18 columns. Collapsing these into one string with
# " ".join(cols[2:]) was the root cause of the false-friend collisions, of
# get_evalue() always returning the KO e-value, and of detecting the AMG
# flag by regex when it is simply column 3. Index them; never join them.
VIB_ANNOT_COLS = [
    "protein", "scaffold", "KO", "AMG", "KO name", "KO evalue", "KO score",
    "KO v-score", "Pfam", "Pfam name", "Pfam evalue", "Pfam score",
    "Pfam v-score", "VOG", "VOG name", "VOG evalue", "VOG score", "VOG v-score",
]

VIB_PROPHAGE_COLS = [
    "scaffold", "fragment", "protein start", "protein stop", "protein length",
    "nucleotide start", "nucleotide stop", "nucleotide length",
]

DF_HMMER_COLS = [
    "hit_id", "replicon", "hit_pos", "hit_sequence_length", "gene_name",
    "i_eval", "hit_score", "hit_profile_cov", "hit_seq_cov",
    "hit_begin_match", "hit_end_match",
]

# ==========================================================================
# 4. OUTPUTS
# ==========================================================================
# vOTU branch
F_VOTU_01 = os.path.join(OUT_DIR, "votu_01_base.tsv")
F_VOTU_02 = os.path.join(OUT_DIR, "votu_02_vibrant.tsv")
F_VOTU_03 = os.path.join(OUT_DIR, "votu_03_coverm.tsv")
F_VOTU_04 = os.path.join(OUT_DIR, "votu_04_phabox_gtdb.tsv")
F_VOTU_05 = os.path.join(OUT_DIR, "votu_05_ct3.tsv")
F_VOTU_06 = os.path.join(OUT_DIR, "votu_06_catalog.tsv")

# gene branch — names kept exactly as the existing scripts already write
# them, so nothing has to be re-run after this config change.
F_GENE_TABLE_L2 = os.path.join(OUT_DIR, "01_gene_table_L2.tsv")
F_GENE_TABLE_L3 = os.path.join(OUT_DIR, "02_gene_table_L3_evidence.tsv")
F_KO_CATS = os.path.join(OUT_DIR, "03_ko_kegg_categories.tsv")
F_GENE_TABLE_CTX = os.path.join(OUT_DIR, "03_gene_table_context.tsv")
F_GENE_TABLE_CLASS = os.path.join(OUT_DIR, "04_gene_table_classified.tsv")

# joined
F_LONG_TABLE = os.path.join(OUT_DIR, "07_long_table_FINAL.tsv")
F_QC_REPORT = os.path.join(OUT_DIR, "99_qc_report.txt")

# ==========================================================================
# 5. SAMPLING DESIGN — the one place fractions are defined
# ==========================================================================
# CB22 / VIRFC : culture-independent (free-living)
# GLI  / LIS   : culture-dependent
# VIRFC / LIS  : MDA-amplified, therefore NOT used for CoverM quantification.
#                Abundance for every vOTU comes from cross-mapping the two
#                non-MDA libraries (CB22, GLI), which is why the TPM columns
#                are TPM_{zone}_{CB22,GLI} only.
# LIS          : mitomycin C induced (prophage induction).
#
# The design is BALANCED for MDA — one MDA and one non-MDA fraction on each
# side — so MDA cannot explain a difference between culture groups.
# It is NOT balanced for mitomycin, which is exclusive to LIS and therefore
# confounded with culture-dependence. CB22 vs GLI (both non-MDA, neither
# induced) is the clean contrast and should be reported as a sensitivity
# analysis of the lifestyle-partitioning result.
FRACTIONS = {
    "CB22":  {"group": "Culture-independent", "mda": "No",  "mitomycin": "No",  "quantified_own": "Yes"},
    "VIRFC": {"group": "Culture-independent", "mda": "Yes", "mitomycin": "No",  "quantified_own": "No"},
    "GLI":   {"group": "Culture-dependent",   "mda": "No",  "mitomycin": "No",  "quantified_own": "Yes"},
    "LIS":   {"group": "Culture-dependent",   "mda": "Yes", "mitomycin": "Yes", "quantified_own": "No"},
}
ZONES = ["DE", "EL", "MB", "OR"]
QUANTIFIED_FRACTIONS = [f for f, v in FRACTIONS.items()
                        if v["quantified_own"] == "Yes"]

# ==========================================================================
# 6. THRESHOLDS
# ==========================================================================
# Same cutoff VIBRANT applies to its own KEGG/Pfam/VOG hmmsearch annotation
# (Kieft, Zhou & Anantharaman 2020, Microbiome 8:90). Citable as-is.
EVALUE_MAX = 1e-5

# Genomic context, DRAM-v style (Shaffer et al. 2020, NAR 48:8883;
# `DRAM-v distill` default is Max auxiliary score: 3). Applied to EVERY ORF,
# so AVGs and defence genes are held to the same criterion — unlike the old
# pipeline, where regex-derived hits were context-filtered and AMG-derived
# ones were not.
CTX_WINDOW_ORFS = 10
CTX_VIRAL_KEYWORDS = [
    "capsid", "tail", "portal", "terminase", "baseplate", "phage",
    "head", "holin", "integrase", "fiber", "virion", "prohead",
]

# Length: a FLAG, never a filter. Pratama et al. 2021 (PeerJ 9:e11447)
# recommend >=10 kb or circular for AMG work; the catalogue keeps 5 kb and
# the conservative subset becomes a filter() in R.
LEN_MAIN_MIN = 5000
LEN_CONSERVATIVE_MIN = 10000

# ==========================================================================
# 7. GENE CLASSIFICATION
# ==========================================================================
# Precedence, applied in genes/05_classify.py:
#   1 specialist_defence : DefenseFinder or PADLOC assigned a defence system
#   2 vibrant_defence    : KO name / Pfam name matches SEARCH_TERMS
#   3 vibrant_amg        : VIBRANT AMG flag (annotations column 3)
#   4 other
# Rule 1 taking precedence over rule 3 is what replaces the old hand-picked
# RECLASSIFY_DEFENSE list of four KOs with a stated criterion.

# Searched ONLY in "KO name" and "Pfam name". Never "VOG name": that field
# describes a reference protein from another organism and was the source of
# every false-friend collision in the old pipeline.
#
# Deliberately dropped from the old SEARCH_TERMS:
#   Innate_Immune_Systems : BREX/DISARM/Gabija/Zorya do not exist as names in
#                           KEGG or Pfam; the old regex returned only dinI
#                           (PF06183). These now come from DefenseFinder and
#                           PADLOC, which found Gabija, Lamassu, Druantia,
#                           Zorya, Pycsar, CBASS, Retron and others.
#   CRISPR_Cas            : the single old hit was uvsW (VOG00935). Real
#                           CRISPR signal (Cas4, 2OG, cas3, cas11, cas12f2)
#                           now comes from PADLOC.
#   bare "acyltransferase": too generic, and redundant with the explicit
#                           "cytolysin-activating" alternative.
SEARCH_TERMS = {
    # AUDITED against the annotations these patterns actually returned. Four
    # terms were removed because they were matching phage STRUCTURAL genes,
    # not virulence factors:
    #
    #   "zona occludens" / "zot"  -- 35 genes. Pfam PF05707 ("Zot") is NOT the
    #     toxin domain: it is the N-terminal ATPase domain of the pI protein,
    #     "the highly conserved protein among filamentous phages", and it is
    #     "essential for the assembly and export of phage virion". Every
    #     filamentous phage carries it. Inoviridae is one of the most abundant
    #     families in this dataset, so this term was counting virion assembly
    #     proteins as host manipulation -- and would have produced enrichment
    #     driven by taxonomic composition rather than by gene function.
    #     Zot's enterotoxic activity lives in the C-terminal domain, which
    #     PF05707 does not cover.
    #
    #   "typhoid toxin" / "ttsA"  -- 25 genes, all carrying Pfam "Glycosyl
    #     hydrolase 108". TtsA "appears to be adapted from endolysins,
    #     phage-related enzymes used by bacteriophages to lyse infected
    #     bacteria and exit from them" (EMBO Rep 2013). On phage contigs a
    #     GH108 hit is the ancestral endolysin, which every phage encodes,
    #     not the Salmonella Typhi toxin-secretion muramidase.
    #
    # What remains are genuine lysogenic-conversion virulence factors.
    # Anchor reference: Waldor & Mekalanos 1996, Science 272:1910,
    # "Lysogenic conversion by a filamentous phage encoding cholera toxin".
    "Toxins_and_Virulence_Factors":
        r"(RTX toxin|cytolysin-activating)",
    "Restriction_Modification":
        r"(mcr[ABC]\b|hsd[RMS]\b|Restriction endonuclease|restriction enzyme)",
    "DNA_Methylation":
        r"(DNA.*methyl|cytosine.*methyl|\bdam\b|\bdcm\b|Hydroxymethylcytosine)",
    "Toxin_Antitoxin_and_Abi":
        r"(Toxin-Antitoxin|TA system|Antitoxin|Maz[EF]\b|Rel[EB]\b|Vap[BC]\b"
        r"|Hok\b|Sok\b|Ccd[AB]\b|Abortive infection|Abi[A-Z]\b|Zeta toxin"
        r"|Yda[ST]\b|YoeB|YefM|Hig[AB]\b|PemK|Gns[AB]\b|SocA)",
    "Superinfection_Exclusion":
        r"(Superinfection exclusion|Sie[AB]\b)",
    # AUDITED: bare "luxR" was matching "LuxR family transcriptional
    # regulator" in any protein carrying that HTH domain, including malT (the
    # maltose regulon activator), which has nothing to do with quorum sensing.
    # The pattern now requires an explicit quorum-sensing context or a named
    # QS regulator. Confirmed by Hampton, Watson & Fineran 2020 (Nature
    # 577:327-336): "Phages also encode LuxR-type proteins, which respond to
    # Gram-negative quorum-sensing signals"; and "These peptide communication
    # systems (which are also known as arbitrium)... inform phages about host
    # availability". See also Silpe & Bassler 2019 (mBio/Cell).
    "Phage_Communication":
        r"(quorum-sensing|bjaR1\b|rpaR\b|rhiR\b|Arbitrium|Aim[RPX]\b)",
    # AUDITED: "tubulin-like" and "cetZ" removed. All 9 hits this category
    # returned were CetZ, an ARCHAEAL tubulin-like cell-shape protein (one hit
    # even carried Pfam "Ribosomal protein S21"). Not one PhuZ or ChmA was
    # found -- and those are the real markers. Phage nucleus formation is
    # restricted to jumbo phages of the family Chimalliviridae, whose genomes
    # run to hundreds of kb; the median contig here is ~7.9 kb, so finding
    # none is the expected result. The category is kept, and returning zero
    # genes is itself a reportable finding.
    #   Nat Struct Mol Biol 2023, s41594-023-01094-5 (ChmA nuclear shell)
    #   Cell Reports 2023, S2211124723004436 (Chimalliviridae core genome)
    "Phage_Nucleus_Assembly":
        r"(PhuZ\b|Chimallin|ChmA\b|Phage nucleus)",
    # REMOVED — "AMR_and_Metal_Resistance" (terB, tehB, beta-lactamase, efflux
    # pump, tellurite resistance). An efflux pump or a beta-lactamase is a
    # metabolic/resistance function the phage carries, not a way of
    # manipulating the host during infection, so calling it "defence" was a
    # category error. These genes now fall through to rule 3 and are counted
    # as AVGs where VIBRANT flagged them as such, which is what they are.
    # ~9 genes; removed on principle, not because of the count.
}

# Kept in the table but NOT counted as anti-phage defence. The reason lives
# here, in the code, so the decision stays visible instead of being lost.
DEFENCE_EXCLUDE = {
    "HEase": "homing endonuclease — a selfish mobile element, ubiquitous in "
             "phage genomes; not an anti-phage defence system",
}

# ==========================================================================
# 7b. DEFENCE CATEGORY VOCABULARY
# ==========================================================================
# The classified table mixes three naming systems in one column: DefenseFinder
# system names (Gabija, Lamassu), PADLOC PROFILE names, which are COMPONENTS
# and not systems (ZorD1 = Zorya component D; DruM2 = Druantia component M),
# raw database accessions used as profile names (cd09637 = Cas4,
# pfam13640 = 2OG), and our own SEARCH_TERMS labels.
#
# The vocabulary below follows the mechanism-based classification of
#   Georjon & Bernheim (2023) Nat Rev Microbiol 21:686-700,
#   "The highly diverse antiphage defence systems of bacteria"
#
# WHY THESE SIX CATEGORIES
# --------------------------------------------------------------------------
# 1. CRISPR-Cas is kept SEPARATE from the other nucleic-acid-degrading
#    systems. It degrades nucleic acid too, but it is the only ADAPTIVE
#    system: it acquires and stores memory of past infections. Lumping it
#    under "degradation" would flatten the one distinction the field is
#    unanimous about.
#
# 2. SIGNALLING IS NOT A SEPARATE CATEGORY. An earlier draft had
#    "signalling-based" alongside "abortive infection". That was wrong.
#    CBASS, Thoeris and Pycsar use cyclic-nucleotide second messengers to
#    SENSE infection, and their effectors then deplete NAD+ and kill the
#    cell. Signalling and abortive infection are one mechanism described at
#    two levels, not two mechanisms.
#
# 3. GENOME MODIFICATION is separated from degradation because the primary
#    act is marking or altering DNA (phosphorothioation by Dnd/Ssp,
#    7-deazapurine by Dpd, methylation), distinguishing self from non-self,
#    rather than cleaving it. BREX and DISARM modify AND degrade; they are
#    assigned to their principal described mechanism and Methods says so,
#    because every review acknowledges the overlap.
#
# 4. HOST MANIPULATION is PHAGE-SIDE and is not bacterial defence at all.
#    Toxins, virulence factors, quorum-sensing regulators, superinfection
#    exclusion and SOS modulation are weapons and controls the phage carries,
#    not immune systems it has captured. Keeping them apart separates "which
#    host immunity does the phage hijack" from "what does the phage do to
#    its host".
#
# 5. "Other / uncharacterised systems" is ONE bucket, not two. An earlier
#    version split "Other characterised systems" from "Uncharacterised /
#    candidate", but after the core-function veto (7c) removed UG39 and MADS
#    what remained on each side was a handful of singletons — Veles, Dodola,
#    MzaA, Menshen, Aristaios, the DS-* series, the PDC-* candidates — none
#    with a mechanism we could report. Two rows in a figure for material we
#    cannot describe is two rows too many. They are counted, named in the
#    supplementary table, and shown as one honest row.
#
# 6. HOMING ENDONUCLEASES (HEase) ARE EXCLUDED FROM DEFENCE ENTIRELY (see
#    DEFENCE_EXCLUDE). They are selfish mobile elements, ubiquitous in phage
#    genomes; counting them as anti-phage defence would inflate the total
#    with something that defends nothing.
# ==========================================================================

# ==========================================================================
# 7b-bis. THE ASSIGNMENT RULE  (added after checking three primary sources)
# ==========================================================================
# Every entry in SYSTEM_TO_CATEGORY below was checked against the full text of:
#
#   Doron et al. 2018,  Science 359:eaar4120        (Table 1: domain content)
#   Hampton, Watson & Fineran 2020, Nature 577:327-336
#   Georjon & Bernheim 2023, Nat Rev Microbiol 21:686-700
#
# THE RULE, stated once so it can go into Methods verbatim:
#
#     A system is assigned to a mechanistic category ONLY where a primary
#     source describes its molecular activity (nuclease, methylase,
#     membrane pore, NAD+ depletion). Systems whose mechanism is not
#     established in the literature are assigned to "Other /
#     uncharacterised systems", regardless of how well characterised their
#     detection is.
#
# This rule can only ever DEFLATE the mechanism categories, never inflate
# them. That is deliberate: the alternative — assigning a system to a
# mechanism on the strength of a plausible-looking domain name — is exactly
# how "DNA polymerase I" ended up being counted as an anti-phage defence
# gene in the previous version of this pipeline.
#
# Georjon & Bernheim report that "the mechanisms of more than half (58%) of
# the reported antiphage systems are not yet deciphered". A large "Other /
# uncharacterised" bucket is therefore the expected outcome of an honest
# classification, not a failure of one. Cite that figure in the Discussion.
#
# APPLYING THE RULE MOVED FOUR SYSTEMS (~98 genes in the catalogue):
#   Druantia  -> Other   Doron Table 1 = "Helicase, methylase". No nuclease.
#   Hachiman  -> Other   Doron Table 1 = "Helicase". No nuclease.
#   Lamassu   -> Other   Doron Table 1 = "SMC ATPase N-terminal domain".
#   Kiwa      -> Other   Doron Table 1 = "No annotated domain".
# None of the six enrichment tests changes sign or loses significance as a
# result; the culture-dependent enrichment survives every category.
#
# STILL RESTING ON A JUDGEMENT CALL, and flagged rather than hidden:
#   - Cas-associated accessory proteins (DEDDh, RT, WYL, PrimPol, 2OG,
#     cd09637) are assigned to CRISPR-Cas on the strength of DefenseFinder's
#     trailing _I_II_III_IV_V_VI type annotation, not a published statement.
#   - BREX and DISARM: both reviews state the mechanism is unresolved. The
#     established part (host-DNA methylation) is what places them.
#   - DNA_Methylation (the largest single label) comes from VIBRANT KO/Pfam
#     name matches with no specialist call, so RM M-subunits and orphan
#     methyltransferases cannot be told apart within it.
#   - Wadjet and Dpd are described as ANTI-PLASMID by Georjon & Bernheim.
# ==========================================================================

DEFENCE_CATEGORIES = [
    "CRISPR-Cas",
    "Nucleic acid degradation",
    "Genome modification",
    "Abortive infection",
    "Other / uncharacterised systems",
    "Host manipulation (phage-encoded)",
]

# Matched against the resolved system name, then the protein name, then the
# raw label; case-insensitively, longest key first.
SYSTEM_TO_CATEGORY = {
    # --- adaptive immunity ------------------------------------------------
    "CRISPR-Cas": "CRISPR-Cas",
    "cas_": "CRISPR-Cas",
    "cas4": "CRISPR-Cas",
    "cas3": "CRISPR-Cas",
    "cas2": "CRISPR-Cas",
    "cas11": "CRISPR-Cas",
    "cas12": "CRISPR-Cas",
    "cas13": "CRISPR-Cas",
    "csf4": "CRISPR-Cas",
    "recD_IV": "CRISPR-Cas",
    # DefenseFinder appends the CRISPR types a gene is found in, so a
    # trailing _I_II_III_IV_V_VI marks a Cas-associated accessory protein.
    "DEDDh_I_II": "CRISPR-Cas",
    "RT_I_II": "CRISPR-Cas",
    "WYL_I_II": "CRISPR-Cas",
    "PrimPol_I_II": "CRISPR-Cas",
    "Specificity_I": "CRISPR-Cas",
    "CARF": "CRISPR-Cas",
    "RtcR": "CRISPR-Cas",
    "2OG": "CRISPR-Cas",
    "mkCas": "CRISPR-Cas",

    # --- innate: cleave the incoming genome -------------------------------
    "RM Type": "Nucleic acid degradation",
    "RM_Type": "Nucleic acid degradation",
    "RM": "Nucleic acid degradation",
    "Restriction-Modification": "Nucleic acid degradation",
    "Restriction_Modification": "Nucleic acid degradation",
    "REase": "Nucleic acid degradation",
    "mREase": "Nucleic acid degradation",
    "HsdR": "Nucleic acid degradation",
    "FS_HsdR": "Nucleic acid degradation",
    "McrBC": "Nucleic acid degradation",
    "VcaM4I": "Nucleic acid degradation",
    "Gabija": "Nucleic acid degradation",      # Doron 2018 Table 1: "ATPase,
                                               # nuclease, helicase"; Georjon &
                                               # Bernheim 2023: "GajA... nonspecific
                                               # DNA nicking endonuclease"
    # Doron 2018 Table 1 lists only "Helicase" for HamAB, because HamA was
    # still a domain of unknown function (DUF1837) in 2018. The mechanism was
    # solved afterwards: HamA IS the effector nuclease -- its DUF1837 domain
    # resembles a PD-(D/E)XK nuclease -- and HamB is the sensor helicase that
    # restrains it until DNA damage is detected. On activation "Hachiman
    # degrades all DNA in the cell".
    #   Cell 2024, "Genome integrity sensing by the broad-spectrum Hachiman
    #     antiphage defense complex" (S0092-8674(24)01068-7)
    #   Nat Commun 2025, "Bacterial Hachiman complex executes DNA cleavage for
    #     antiphage defense" (s41467-025-57851-1)
    "Hachiman": "Nucleic acid degradation",
    "Shedu": "Nucleic acid degradation",       # Doron 2018 Table 1: "Nuclease"
    "Wadjet": "Nucleic acid degradation",      # Doron 2018: "MukBEF condensin,
                                               # topoisomerase VI". NOTE: Georjon &
                                               # Bernheim call Wadjet ANTI-PLASMID.
    # MOVED to "Other" — Doron 2018 Table 1 lists "Helicase, methylase" for
    # Druantia, plus DUF1998 and "highly variable genes with no recognizable
    # domains or function prediction". There is NO nuclease domain anywhere in
    # the primary description. This is the single largest reassignment made
    # under the rule above (~225 genes) and it DEFLATES the headline
    # "Nucleic acid degradation" count, which is the honest direction.
    "Druantia": "Other / uncharacterised systems",
    "Septu": "Nucleic acid degradation",       # Doron 2018 Table 1: "ATPase,
                                               # HNH nuclease"
    "Ptu": "Nucleic acid degradation",
    "Aga_Nuclease": "Nucleic acid degradation",
    "Nuclease": "Nucleic acid degradation",

    # --- innate: mark or alter the DNA ------------------------------------
    "DNA methylation": "Genome modification",
    "DNA_Methylation": "Genome modification",
    "MTase": "Genome modification",
    "oMTase": "Genome modification",
    "DNA-modification": "Genome modification",
    "Dpd": "Genome modification",              # Georjon & Bernheim: "inserts
                                               # 7-deazaguanine derivatives in the
                                               # host DNA to distinguish it from the
                                               # phage invading and non-modified
                                               # DNA". Same review labels it
                                               # anti-plasmid; note in Methods.
                                               # SEE 7c: the diagnostic gene, dpdA,
                                               # is vetoed OUT of this list on
                                               # purpose -- it is real evidence.
                                               # folE/queC/queD/queE are vetoed
                                               # IN 7c because DefenseFinder's own
                                               # Dpd model needs dpdA-dpdK besides
                                               # these four, and none of the "Dpd"
                                               # calls in this dataset ever include
                                               # it (see 7c for the numbers).
    "Dnd": "Genome modification",              # G&B: "DndACDE... double-stranded
                                               # phosphorothioate modification of
                                               # the host DNA"
    "Ssp": "Genome modification",              # G&B: "SspABCD... single-stranded
                                               # phosphorothioate modification"
    "BREX": "Genome modification",
    "DISARM": "Genome modification",
    "Drm": "Genome modification",
    "PhiR12.1_MorcDcm": "Genome modification",

    # --- innate: kill the cell to save the population ---------------------
    "CBASS": "Abortive infection",
    "Thoeris": "Abortive infection",
    "Ths": "Abortive infection",
    "Pycsar": "Abortive infection",
    "Pyc_": "Abortive infection",
    "Toxin-Antitoxin": "Abortive infection",
    "Toxin_Antitoxin": "Abortive infection",
    "MazEF": "Abortive infection",
    "HicAB": "Abortive infection",
    "mEp506": "Abortive infection",
    "PhiR47.1": "Abortive infection",
    "RosmerTA": "Abortive infection",
    "PsyrTA": "Abortive infection",
    "ShosTA": "Abortive infection",
    "DarTG": "Abortive infection",
    "PanTA": "Abortive infection",
    "Abi": "Abortive infection",
    # Doron 2018 gave only "SMC ATPase N-terminal domain". The structure and
    # activation mechanism were solved in 2025: LmuB (SMC family) is the
    # sensor, and LmuA -- described by the authors as "the Abi-inducing
    # effector" -- is released on activation, tetramerises, and cleaves DNA,
    # "halting phage propagation by terminating the infected host".
    #   Nat Struct Mol Biol 2025, "Structure and activation mechanism of a
    #     Lamassu phage and plasmid defense system" (s41594-025-01677-4)
    # AMBIGUITY, stated rather than hidden: the effector is a nuclease, so an
    # effector-only reading would put Lamassu under nucleic acid degradation.
    # It sits under Abortive infection because the source names the OUTCOME
    # (Abi, host death) explicitly, and because LmuA effectors are
    # enzymatically diverse (nuclease, protease) while the Abi outcome is
    # constant. Gabija is the mirror case and is filed by its nuclease.
    "Lamassu": "Abortive infection",
    "Zorya": "Abortive infection",             # TRIPLE confirmed: Doron 2018
                                               # ("MotA/ExbB, MotB", "Membrane
                                               # associated"); Hampton 2020
                                               # ("proposed to be an Abi system...
                                               # depolarization and cell death");
                                               # Georjon & Bernheim 2023
    # MOVED to "Other" — Doron 2018 Table 1: "No annotated domain". The
    # "Membrane associated" comment hints at an Abi mechanism but no source
    # states one.
    "Kiwa": "Other / uncharacterised systems",
    "Retron": "Abortive infection",
    "retron": "Abortive infection",
    "DRT": "Abortive infection",
    "drt": "Abortive infection",
    "SIR2": "Abortive infection",
    "TIR": "Abortive infection",
    "SEFIR": "Abortive infection",
    "dCTPdeaminase": "Abortive infection",     # Tal et al. 2022, Nat Microbiol
                                               # s41564-022-01158-0: these enzymes
                                               # "completely eliminate the specific
                                               # deoxynucleotide... thus starving
                                               # the phage of an essential DNA
                                               # building block and halting its
                                               # replication". NOTE: no host death
                                               # is described, so this is not Abi in
                                               # the strict sense. Kept here only
                                               # because the vocabulary has no
                                               # nucleotide-depletion category and
                                               # n = 1. NOTE: Georjon & Bernheim give
                                               # nucleotide-pool depletion its own
                                               # effector class rather than filing
                                               # it under Abi. Kept here because
                                               # depletion arrests/kills the cell;
                                               # state the departure in Methods.
    "SoFic": "Abortive infection",
    "IetA": "Abortive infection",
    "Gao_Iet": "Abortive infection",
    "QatA": "Abortive infection",
    "Gao_Qat": "Abortive infection",

    # --- described or candidate, no reportable mechanism ------------------
    # --- described mechanism, therefore NOT "Other" -----------------------
    # Georjon & Bernheim 2023: "Bacterial gasdermins (bGSDMs) are generally
    # found next to one or several proteases, which are activated during phage
    # infection to cleave an inhibitory C-terminal peptide of the bGSDMs...
    # Activated bGSDMs assemble and form large pores in the membrane,
    # eventually killing the host."  Membrane pore -> cell death -> Abi.
    #
    # NOTE ON WHY THIS KEY IS NEEDED AT ALL: PADLOC's raw profile name for
    # these genes is "subtilase" (the protease partner), and "subtilase" is
    # listed below under "Other". to_category() tries the RESOLVED SYSTEM NAME
    # before the raw label, so "Gasdermin" is now matched first and the
    # "subtilase" entry never fires for them. Order within this dict is
    # irrelevant here -- it is the candidate order in to_category() that
    # decides. The "subtilase" key stays as a fallback for any gene whose
    # PADLOC yaml lookup fails.
    "Gasdermin": "Abortive infection",

    "Dodola": "Other / uncharacterised systems",
    "Veles": "Other / uncharacterised systems",
    "Menshen": "Other / uncharacterised systems",
    "Aristaios": "Other / uncharacterised systems",
    "Sirona": "Other / uncharacterised systems",
    "Taranis": "Other / uncharacterised systems",
    "Sucellos": "Other / uncharacterised systems",
    "Reve": "Other / uncharacterised systems",
    "Dag": "Other / uncharacterised systems",
    "BstA": "Other / uncharacterised systems",
    "RAZR": "Other / uncharacterised systems",
    "TagI": "Other / uncharacterised systems",
    "Shango": "Other / uncharacterised systems",
    "Hailong": "Other / uncharacterised systems",
    "Paris": "Other / uncharacterised systems",
    "Brig": "Other / uncharacterised systems",
    "Aditi": "Other / uncharacterised systems",
    "Old": "Other / uncharacterised systems",
    # MADS was uncharacterised when this vocabulary was first written. Cell
    # Host & Microbe 2024 describes it as "eight genes with predicted nuclease,
    # ATPase, kinase, and methyltransferase domains, most of which are
    # essential for either self/non-self discrimination, DNA restriction, or
    # both", and notes it is "predicted to encode an N6 methyladenosine MTase,
    # which suggests that self/non-self discrimination relies on DNA
    # methylation". That is the BREX/DISARM pattern, so it is filed with them.
    #   Cell Host Microbe 2024, "The bacterial defense system MADS interacts
    #     with CRISPR-Cas to limit phage infection and escape"
    #     (S1931-3128(24)00264-6)
    # JUDGEMENT CALL, flagged: the same source lists a nuclease domain, so
    # "Nucleic acid degradation" is defensible too. Genome modification is
    # chosen for consistency with BREX and DISARM, which are placed by their
    # methylation-based discrimination rather than by their nuclease half.
    # Note that PADLOC's own yaml already groups MADS under "DNA-modification
    # systems (other)"; this restores that call.
    "MADS": "Genome modification",
    "mads": "Genome modification",
    "Dynamins": "Other / uncharacterised systems",
    "DinG": "Other / uncharacterised systems",
    "PbeA": "Other / uncharacterised systems",
    "MzaA": "Other / uncharacterised systems",
    "MzaB": "Other / uncharacterised systems",
    "Mza": "Other / uncharacterised systems",
    "SNIPE": "Other / uncharacterised systems",
    "Effector": "Other / uncharacterised systems",
    "subtilase": "Other / uncharacterised systems",
    "Butters_gp30": "Other / uncharacterised systems",
    "MMB_gp29": "Other / uncharacterised systems",
    "Armada": "Other / uncharacterised systems",
    "VP1839": "Other / uncharacterised systems",
    "SDIC": "Other / uncharacterised systems",
    "UG": "Other / uncharacterised systems",
    "DS-": "Other / uncharacterised systems",
    "PDC": "Other / uncharacterised systems",
    "PD-": "Other / uncharacterised systems",
    "HEC": "Other / uncharacterised systems",
    "cluster": "Other / uncharacterised systems",
    "icity": "Other / uncharacterised systems",
    "COG": "Other / uncharacterised systems",
    "cd0": "Other / uncharacterised systems",
    "pfam": "Other / uncharacterised systems",
    "FS_": "Other / uncharacterised systems",
    "HTH": "Other / uncharacterised systems",
    "E1E2": "Other / uncharacterised systems",
    "V_U": "Other / uncharacterised systems",
    "Type_V": "Other / uncharacterised systems",
    "Control_protein": "Other / uncharacterised systems",

    # --- phage-side, NOT bacterial defence --------------------------------
    "Toxins and virulence": "Host manipulation (phage-encoded)",
    "Toxins_and_Virulence_Factors": "Host manipulation (phage-encoded)",
    "Phage communication": "Host manipulation (phage-encoded)",
    "Phage_Communication": "Host manipulation (phage-encoded)",
    "Superinfection exclusion": "Host manipulation (phage-encoded)",
    "Superinfection_Exclusion": "Host manipulation (phage-encoded)",
    "Phage nucleus": "Host manipulation (phage-encoded)",
    "Phage_Nucleus_Assembly": "Host manipulation (phage-encoded)",
    # REMOVED — dead key. No SEARCH_TERM could ever produce the label
    # "SOS_response_modulation", so this entry was unreachable. dinI
    # (PF06183) is a HOST SOS gene (81 aa, inhibits the RecA coprotease
    # so LexA is not cleaved), not a phage strategy; the characterised
    # phage-encoded SOS modulator is a different protein, gp7 of phage
    # GIL01 (NAR 2015). A dinI homologue on a phage contig is more
    # likely host-derived carry-over, so it is deliberately NOT counted.
}


# ==========================================================================
# 7b-ter. LABEL_OVERRIDE — raw labels whose automatic resolution is wrong
# ==========================================================================
# 06_harmonise_categories.py resolves a raw label by trying, in order, the
# system name from PADLOC's yaml, then the protein names from hmm_meta, then
# the raw label itself. That order is right almost always, and it is what
# rescues meaningless accessions like cd09637 by looking up what protein the
# profile actually describes.
#
# It is wrong for one specific family of labels: DefenseFinder and PADLOC
# annotate DOMAINS that CO-OCCUR with CRISPR-Cas systems, and name them after
# the CRISPR types they are seen in. "DEDDh_I_II_III_IV_V_VI" means "a DEDDh
# exonuclease domain observed in CRISPR types I to VI". It does not mean the
# protein is a Cas protein, and there is no such thing as a "DEDDh system".
#
# The audit made this unambiguous:
#   * 0 of 176 genes in the CRISPR-Cas category belonged to a COMPLETE
#     system (df_in_complete_system / pl_in_complete_system both 0).
#   * Only 2 of 176 carried an annotation naming a Cas protein.
#   * The bulk were generic exonucleases (MGME1, RNase T, RNase H), group II
#     intron reverse transcriptases (ltrA), and Tn7/Mu transposition proteins
#     (TnsC, MuB) -- none of which is CRISPR-Cas by any reading.
#
# This config already anticipated it without connecting the dots, in 7c:
#   "With only 1 of 166 cd09637 hits actually annotated cas4, erring towards
#    rejection is the right trade."
#
# Three of these labels cannot be fixed by editing SYSTEM_TO_CATEGORY,
# because they never reach the raw-label step: cd09637 and COG1468 resolve
# through hmm_meta's protein.name ("Cas4"), and TnsC through PADLOC's yaml
# ("CRISPR-Cas Type I-F3"). Hence an explicit override, consulted first.
#
# THE RULE APPLIED HERE, same as everywhere else in this file: a gene is
# assigned to a mechanistic category only when it belongs to a system with a
# described mechanism. A domain seen near CRISPR systems is not a system.
# Labels naming an actual Cas gene (cas2, cas3, cas4, cas11, cas12, csf4)
# are deliberately NOT overridden and stay in CRISPR-Cas.
#
# Matched by longest prefix, case-insensitively, against the RAW label only.
# ==========================================================================
LABEL_OVERRIDE = {
    # --- exonuclease / RT / helicase domains named after CRISPR types -----
    "DEDDh": "Other / uncharacterised systems",
    "RT_I_II": "Other / uncharacterised systems",
    "PrimPol": "Other / uncharacterised systems",
    "WYL": "Other / uncharacterised systems",
    "recD_IV": "Other / uncharacterised systems",
    "DinG": "Other / uncharacterised systems",
    # --- database accessions that hmm_meta maps to "Cas4" -----------------
    # config 7c already records that only 1 of 166 cd09637 hits was really
    # annotated cas4; the rest are RecB-family exonucleases.
    "cd09637": "Other / uncharacterised systems",
    "COG1468": "Other / uncharacterised systems",
    "pfam13640": "Other / uncharacterised systems",
    "pfam01930": "Other / uncharacterised systems",
    # --- transposition machinery, not immunity ----------------------------
    "TnsC": "Other / uncharacterised systems",
    # --- CARF/RtcR signalling domains and unnamed accessions --------------
    "RtcR": "Other / uncharacterised systems",
    "CARF": "Other / uncharacterised systems",
    "shah_acc": "Other / uncharacterised systems",
    "icity": "Other / uncharacterised systems",
    "mkCas": "Other / uncharacterised systems",
}

# ==========================================================================
# 7c. CORE-FUNCTION VETO
# ==========================================================================
# Cross-referencing every specialist (PADLOC / DefenseFinder) call against
# the protein's OWN VIBRANT annotation showed that a quarter of them were
# core cellular machinery captured by the HMM profiles of poorly
# characterised defence families. Evidence, from the classified gene table:
#
#   label       n     what the protein's own KO/Pfam actually says
#   ---------------------------------------------------------------------
#   UG9        609    527x Pfam "DNA polymerase family A", 523x KO polA
#   UG6        166    Toprim + DnaB-like helicase; KO dnaG, PEO1 (primase)
#   UG39        75    69x Pfam "VirE N-terminal domain" = primase
#   DEDDh_*    238    75x KO dnaQ = DNA Pol III epsilon (proofreading)
#   cd09637    166    98x KO K07465 RecB exonuclease; only 1 was cas4
#   pfam13640  163    110x KO K07336 / P4HA prolyl hydroxylase (metabolic)
#   SspBCDE    168    37x KO cysH (PAPS reductase) vs only 3x dndC
#   Dodola     129    AAA+ ATPases: clpA, clpC, cobS
#   MADS        73    17x KO alpA (prophage regulator), 17x HTH, median 81 aa
#   RM          96    58x radD, 24x IRC3, 83x VOG uvsW
#   ZorD1      458    SNF2 helicases (SMARCAL1); Lamassu 231 -> SbcC/Rad50
#   Dpd         83    "folE"/"queC"/"queD"/"queE" -- the shared preQ0/queuosine
#               biosynthesis genes, in EVERY case. DefenseFinder's own Dpd model
#               requires 15 genes (FolE, QueD, QueE, QueC, plus DpdA through
#               DpdK -- see defensefinder.mdmlab.fr/wiki/defense-systems/dpd).
#               DpdA is the transglycosylase that actually inserts the modified
#               base into DNA; without it there is no evidence the DNA is
#               modified at all, only that the housekeeping preQ0 pathway is
#               present -- which VIBRANT already calls correctly as an AVG
#               (Huang et al. 2021, Environ Microbiol, PMID 33511765, lists the
#               same four genes as roseophage AMGs). In this dataset: 0 of 83
#               Dpd-tagged genes belong to a complete system
#               (df_in_complete_system / pl_in_complete_system both 0 for all
#               83), and the largest per-genome cluster is 4 genes -- exactly
#               folE+queC+queD+queE, never dpdA. Same failure mode as the
#               CRISPR-Cas type-suffix problem in 7b-ter: a system name
#               assigned on partial, non-diagnostic evidence.
#
# The RM row is the SAME false-friend problem the original pipeline had:
# radD, IRC3 and uvsW are repair helicases sharing PF04851 "Type III
# restriction enzyme, res subunit". The original FALSOS_AMIGOS blacklist
# excluded them correctly; PADLOC's generic "RM" profile pulls them back in.
#
# THE RULE
#   A specialist defence assignment is REJECTED when the protein's own
#   annotation identifies it as a core replication, repair or metabolic
#   enzyme, or as a prophage lifecycle regulator.
#
# BOTH IDENTITY FIELDS ARE SEARCHED. An earlier version read only 'KO name',
# which let UG39 through: 69 of its 75 proteins carry the VirE primase
# domain but have no KO assignment at all. Rule 2 already searched KO name
# AND Pfam name; the veto now does the same, which is the consistent
# behaviour.
#
# WHY IT IS APPLIED TO ALL SPECIALIST CALLS, not only to uncharacterised
# families: simpler to state, needs no list of which families count as
# "characterised", and errs towards rejection. In practice named systems
# (Gabija, Zorya, CBASS) are almost never hit, because their profiles are
# specific enough that their proteins do not carry a polA or dnaQ orthologue.
#
# THE KNOWN COST, stated rather than hidden: Cas4 IS a RecB-family
# exonuclease, so rejecting K07465 discards a handful of genuine Cas4 calls
# along with ~98 generic exonucleases. With only 1 of 166 cd09637 hits
# actually annotated cas4, erring towards rejection is the right trade.
#
# NOTHING IS DISCARDED. A vetoed gene falls through to rules 2 and 3, so a
# cysH or cobS rejected here is picked up as an AVG, which is what it is.
# Every veto is counted and written into class_source and specialist_vetoed.
#
# Methods sentence:
#   "Defence assignments made by PADLOC or DefenseFinder were rejected when
#    the protein's own KEGG orthologue or Pfam domain identified it as a
#    core replication, repair or metabolic enzyme, or as a prophage
#    lifecycle regulator (Supplementary Table SX); such genes were
#    reassigned to the auxiliary viral gene set where applicable."
# ==========================================================================

CORE_FUNCTION_REJECT = (
    r"("
    # --- replication ----------------------------------------------------
    r"\bpolA\b|DNA polymerase I\b|DNA polymerase family A"
    r"|\bdnaQ\b|DNA polymerase III"
    r"|\bdnaG\b|DNA primase"
    r"|\bdnaB\b|replicative DNA helicase|DnaB-like helicase"
    r"|\bPEO1\b|twinkle"
    # UG39: 69 of 75 proteins are "VirE N-terminal domain", a primase
    # domain, and carry no KO at all — only reachable via Pfam.
    r"|VirE N-terminal|VirE_N|Toprim"
    # --- repair helicases: the original false friends ---------------------
    r"|\bradD\b|\bIRC3\b|\bradA\b|\bsms\b"
    r"|\bsbcC\b|\bsbcD\b|\brad50\b|\bmre11\b"
    r"|\bERCC3\b|\bERCC4\b|\bSMARCAL1\b|SNF2 family"
    # --- generic nucleases that are not, by themselves, defence -----------
    r"|RecB family exonuclease|K07465"
    # Cas4 IS a PD-(D/E)XK / RecB-family exonuclease, so the Cas4 profiles hit
    # every other exonuclease of that fold. Listing the raw labels confirmed it:
    # of the 25 genes in CRISPR-Cas, 13 were MGME1 (mitochondrial genome
    # maintenance exonuclease), 3 addB, 1 dinG, 1 cbf -- all under labels
    # literally named Cas4_II / cas4_I_II_V_maka / csf4_IV-A / cas3_HD. Only 2
    # carried the annotation "cas4; CRISPR-associated exonuclease Cas4".
    #
    # THIS CANNOT BE FIXED BY LABEL_OVERRIDE: the two genuine cas4 genes carry
    # the SAME raw label as six of the false ones. The protein's own annotation
    # is the only thing that separates them, which is precisely what this veto
    # is for -- the same rule that already rejects radD and IRC3 under an "RM"
    # label. MGME1 and AddB are genome-maintenance/repair nucleases, DinG is a
    # repair helicase: all core function by the rule stated above.
    r"|\bMGME1\b|\bDDK1\b|genome maintenance exonuclease"
    r"|\baddB\b|helicase/nuclease subunit B"
    r"|\bdinG\b"
    r"|\bcbf1?\b|3'-5' exoribonuclease"
    # Cas12f (Cas14) is derived from IS200/IS605 TnpB transposases, so the
    # cas12f2 profile hits transposases. Transposition machinery is not
    # immunity; TnsC and MuB are already handled the same way in 7b-ter.
    r"|K07496|utative transposase|robable transposase"
    # --- DNA repair by dealkylation: a repair enzyme, not a methylase.
    # AlkB was slipping into "Genome modification" because SEARCH_TERMS'
    # DNA_Methylation pattern "DNA.*methyl" matches "DNA oxidative
    # DEmethylase". Only 1 gene here, added for consistency with radD/IRC3.
    r"|\balkB\b|oxidative demethylase"
    # --- preQ0/queuosine biosynthesis: shared with the universal AMG pathway,
    # not diagnostic of Dpd on its own. dpdA is deliberately NOT here -- see
    # the Dpd entry above and the evidence table for the 83-gene count.
    r"|\bfolE\b|GTP cyclohydrolase"
    r"|\bqueC\b|7-cyano-7-deazaguanine synthase|Queuosine biosynthesis protein QueC"
    r"|\bqueD\b|6-pyruvoyl.{0,4}tetrahydropterin"
    r"|\bqueE\b|7-carboxy-7-deazaguanine synthase"
    # --- prophage lifecycle regulation is not immunity --------------------
    # MADS: 17x KO alpA, 17x HTH domain, median protein length 81 aa.
    # AlpA regulates the prophage's own excision.
    r"|\balpA\b|prophage regulatory|Prophage CP4-57"
    # --- metabolism (several are VIBRANT AMGs in their own right) ---------
    r"|PKHD-type hydroxylase|prolyl 4-hydroxylase|K07336|\bP4HA\b"
    r"|\bcysH\b|phosphoadenosine phosphosulfate reductase"
    r"|\bclpA\b|\bclpC\b|ATP-dependent Clp protease"
    r"|\bcobS\b|cobaltochelatase"
    r")"
)

# ==========================================================================
# 8. HELPERS
# ==========================================================================
def check_header(path, expected, sep="\t", label=None):
    """Assert the file's header matches the expected schema.

    Positional column indices are only safe if the header is verified. A
    silently shifted column yields plausible wrong numbers, which is the
    worst possible failure mode.
    """
    label = label or os.path.basename(path)
    with open(path) as fh:
        got = [c.strip() for c in fh.readline().rstrip("\n").split(sep)]
    if got[:len(expected)] != expected:
        raise SystemExit(
            f"[SCHEMA ERROR] {label}\n"
            f"  expected: {expected}\n"
            f"  found   : {got[:len(expected)]}")
    return got


def require(*paths):
    """Fail immediately, and loudly, if an input file is missing."""
    missing = [p for p in paths if not os.path.exists(p)]
    if missing:
        raise SystemExit("[MISSING INPUT]\n  " + "\n  ".join(missing))


class Counter_:
    """Running tally printed at the end of every script.

    Every exclusion, unmatched merge and fallback increments a named
    counter, so nothing is ever lost in silence.
    """

    def __init__(self, title):
        self.title = title
        self.d = {}

    def add(self, key, n=1):
        self.d[key] = self.d.get(key, 0) + n

    def report(self, total=None):
        print("\n" + "-" * 66)
        print(f" {self.title}")
        print("-" * 66)
        if not self.d:
            print("   (nothing to report — clean run)")
        for k, v in sorted(self.d.items(), key=lambda kv: -kv[1]):
            pct = f"  ({100.0 * v / total:.1f}%)" if total else ""
            print(f"   {v:>9}  {k}{pct}")
