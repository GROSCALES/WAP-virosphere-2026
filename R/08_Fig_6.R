# =============================================================================
# 08_Fig6.R - Figure 6. Two temperate genomes, opposite fractions, opposite
# arsenals. THE CONTROL FIGURE.
#
# CHANGES 23/9
#   * The title of panel b gives the vOTU length from Data S1 (42.8 kb), not
#     the length of the contig (43.2 kb). Panel a was already correct.
#   * MAPS_DIR: the absolute path lives in one variable now.
#   * Panel b is titled Peduoviridae, the family Data S1 assigns to that vOTU;
#     panel a has no family assignment and keeps the class.
#
# WHY THIS FIGURE EXISTS
#   Figures 4 and 5 change two variables at once: Figure 4 is virulent and
#   culture-independent, Figure 5 is temperate and culture-dependent. A reader
#   comparing them cannot tell whether the difference in arsenal follows the
#   fraction or the lifestyle -- and the statistics say it is the fraction
#   (lifestyle held constant: CMH OR = 1.01, p = 0.98; fraction held constant:
#   OR = 3.89, p = 8.9e-11). This figure holds the lifestyle constant instead:
#   both genomes are Temperate (predicted), both are 100% complete by CheckV,
#   and only the fraction differs.
#
#        panel a   VIRFC_NODE_3   VIR, culture-independent   2/2 self-protection
#        panel b   EL-GLI_NODE_119  CUL, culture-dependent   2/2 exclusion
#
# GENES ARE MATCHED BY COORDINATE, NOT BY ORF NUMBER
#   The .gbf files come from Cenote-Taker 3, which calls genes with
#   pyrodigal-gv; the pipeline's orf_num comes from Prodigal 2.6.3 run with
#   -p meta. The two numberings need not agree, so the gene assignments below
#   carry the START/END from the pipeline tables and are matched to the parsed
#   GenBank record by overlap. Section 5 prints every match and every miss:
#   check that output before trusting the figure.
#
# WHAT IS LABELLED, AND WHAT IS NOT
#   Panel a: the four methyltransferases, plus the cysH/SspD pair. The three
#   tandem methylases (17,294-21,411, contiguous) were each called an RM Type
#   II component by DefenseFinder AND PADLOC; there is no restriction
#   endonuclease anywhere in the 65.6 kb contig, so the reconstruction
#   reassigns them as orphan methyltransferases. They are coloured as Genome
#   modification, which is where the reconstruction puts them, not as Nucleic
#   acid degradation, which is the raw gene-level call. That reassignment is
#   the point of the panel and it is stated in the legend.
#
#   Panel b: gajA and gajB only. Two further defence genes are drawn in their
#   category colour but deliberately left unlabelled:
#     - the dam at 15,073-15,924, which PADLOC itself calls oMTase. The
#       reconstruction scores an RM Type II locus here by pairing it with GajA
#       8 kb away; that is one of the nine cases Supplementary Methods S11
#       flags as the least secure calls in the whole reconstruction, so the
#       figure makes no claim about it.
#     - the CI-repressor-domain protein at 12,508-13,140, a DS-34 candidate
#       whose locus never reconstructed an architecture.
#   Same rule as luxR and the transposases in Figure 5: colour shows the
#   category, a label would assert a system.
#
# THE BRACKET IS ONLY DRAWN OVER A COMPLETE, CONFIRMED LOCUS
#   L00454 (Gabija, GajA + GajB, 2/2 components, confidence high) is the ONLY
#   complete Gabija among the 202 reconstructed in the whole catalogue. The
#   other 201 fail because GajB is by definition a UvrD-family helicase and a
#   lone GajB call cannot be told from housekeeping. Only the pair is
#   informative, and it happens once.
# =============================================================================

library(ggplot2)
library(dplyr)
library(gggenes)
library(cowplot)

# ------------------------------------------------------------------------------
# 1. PATHS AND METADATA
# ------------------------------------------------------------------------------
# One line to change for the repository: MAPS_DIR holds the two .gbf files and
# palette_maps.R, and it is the only absolute path in the script.
MAPS_DIR <- "C:/Users/gasch/Documents/R/vOTUs-R-V2"

path_a <- file.path(MAPS_DIR, "vOTU_viral_51405.gbf")   # MCHA-DE-14-VIRFC_NODE_3
path_b <- file.path(MAPS_DIR, "vOTU_viral_7328.gbf")    # MCHA-EL-19-GLI_NODE_119

NAME_a <- "MCHA-DE-14-VIRFC_NODE_3"
NAME_b <- "MCHA-EL-19-GLI_NODE_119"

# setwd(MAPS_DIR) USED TO BE HERE, and it is the one thing outside the export
# block that changed. It did two harmful things and one necessary one:
#   harmful  - it left the whole R session pointing at MAPS_DIR, so every
#              script run afterwards looked for its tables, its palette and
#              its output folder in the wrong place;
#   harmful  - it sent this figure's files to MAPS_DIR instead of your R
#              folder;
#   necessary - it made source("palette_maps.R") below resolve to the copy
#              in MAPS_DIR, which is THE PALETTE THIS FIGURE WAS DRAWN WITH.
# The line below keeps the necessary part and drops the harmful parts: the
# same palette file is read, by its full path, and the working directory is
# left alone. The figure is pixel-identical to the one you approved.
source(file.path(MAPS_DIR, "palette_maps.R"))
paleta_global <- MAP_PALETTE

SZ_GENE  <- 5.0
SZ_AXIS  <- 13
SZ_TITLE <- 16
SZ_NOTE  <- 4.3
# The floor used to be 0.88 and the bracket caption sat at 0.885, so the text
# was drawn straight onto the axis line. There is now room below the arrows
# for the bracket and its caption.
Y_MAP    <- c(0.66, 1.52)

# THE dam AT 15,073 IS NOW LABELLED. It was left blank because the
# reconstruction pairs it with GajA 8 kb away to score an RM Type II, and that
# is one of the nine calls Supplementary Methods S11 flags as least secure --
# so the figure must not assert that system. But `dam` is a GENE name, not a
# system name: writing it asserts "adenine methyltransferase", which is what
# PADLOC itself calls it (oMTase), and it is exactly the standard already
# applied to dcm, sspD and the N6/N4-MTases in panel a. What must not happen
# is a BRACKET spanning it and gajA, and there is none. An unlabelled purple
# gene 6 kb from two labelled purple genes reads as a third Gabija component,
# which is worse than the risk it was avoiding. Set to FALSE to restore the
# blank; nothing else changes.
LABEL_DAM_B <- TRUE

# Labels at 45 degrees collided wherever two annotated genes were adjacent --
# cysH/sspD overlap by 4 bp, and the three methyltransferases are contiguous.
# Same generic solution as Figure 5: horizontal labels stacked on the lowest
# free level, joined to their gene by a leader line. Nothing is hand-placed,
# so it cannot break if a coordinate changes.
asignar_niveles <- function(df, min_dist, niveles = c(1.16, 1.28, 1.40)) {
  df$label_y <- NA_real_
  idx <- which(df$label != "")
  if (!length(idx)) return(df)
  idx <- idx[order(df$start[idx])]
  ultimo <- rep(-Inf, length(niveles))
  for (i in idx) {
    libres <- which(df$start[i] - ultimo >= min_dist)
    k <- if (length(libres)) min(libres) else which.min(ultimo)
    df$label_y[i] <- niveles[k]
    ultimo[k] <- df$start[i]
  }
  df
}

tema_isme <- theme_minimal() + theme(
  panel.grid = element_blank(), axis.line.x = element_line(color = "black"),
  axis.text.y = element_blank(), axis.title.y = element_blank(),
  axis.title.x = element_text(size = SZ_AXIS - 1),
  axis.text.x = element_text(size = SZ_AXIS, colour = "grey20"),
  axis.ticks.x = element_line(colour = "grey20", linewidth = 0.4),
  legend.position = "none",
  plot.title = element_text(face = "bold", size = SZ_TITLE),
  plot.margin = margin(t = 12, r = 20, b = 16, l = 20)
)

# Panel b carries a family-level assignment (Peduoviridae) and panel a does
# not, so each title names the deepest rank that vOTU actually has.
ttl6 <- function(letter, name, kb, zone, frac, origin, taxon = "Caudoviricetes") {
  bquote(bold(.(paste0(letter, "  Temperate "))) *
           bolditalic(.(taxon)) *
           bold(.(paste0(" \u00b7 ", name, " (", kb, " \u00b7 ", zone,
                         " \u00b7 ", frac, " \u00b7 ", origin, ")"))))
}

# ------------------------------------------------------------------------------
# 2. PARSER  (same as 07_Fig4.R)
# ------------------------------------------------------------------------------
parse_gbk <- function(file) {
  lines <- readLines(file)
  cds_idx <- grep("^     CDS ", lines)
  data <- do.call(rbind, lapply(cds_idx, function(i) {
    loc_line <- lines[i]
    nums <- as.numeric(regmatches(loc_line, gregexpr("[0-9]+", loc_line))[[1]])
    prod <- "hypothetical protein"
    for (j in (i + 1):min(i + 20, length(lines))) {
      if (grepl("/product=", lines[j])) {
        prod <- gsub('.*/product="([^"]+)".*', "\\1", lines[j]); break
      }
    }
    data.frame(start = min(nums), end = max(nums),
               strand = ifelse(grepl("complement", loc_line), -1, 1),
               product = prod)
  }))
  data$gene_num <- seq_len(nrow(data))
  data$direccion_hebra <- data$strand == 1
  data
}

# ------------------------------------------------------------------------------
# 3. GENE ASSIGNMENTS - coordinates taken from 07_long_table_FINAL_28_8.tsv
#    Category names are the ones palette_maps.R already understands.
# ------------------------------------------------------------------------------
ann_a <- tribble(
  ~start,  ~end,   ~label,      ~cat,
  1555,   2325,  "dcm",        "Defence: Genome modification",
  17294,  18625,  "N6-MTase",   "Defence: Genome modification",   # see header
  18668,  20419,  "dcm",        "Defence: Genome modification",
  20419,  21411,  "N4-MTase",   "Defence: Genome modification",
  25773,  26780,  "cysH",       "AMG: Energy metabolism",
  26777,  27577,  "sspD",       "Defence: Genome modification"
)
ann_b <- tribble(
  ~start,  ~end,   ~label,      ~cat,
  4947,   7106,  "gajA",       "Defence: Nucleic acid degradation",
  7109,   8860,  "gajB",       "Defence: Nucleic acid degradation",
  12508,  13140,  "",           "Defence: Other / uncharacterised",  # unlabelled on purpose
  15073,  15924,  "dam",        "Defence: Nucleic acid degradation"  # see LABEL_DAM_B
)
if (!LABEL_DAM_B) ann_b$label[ann_b$start == 15073] <- ""

# cysH keeps its real category. palette_maps.R now carries
# "AMG: Energy metabolism"; it used to be recoloured here as cofactors and
# vitamins, which made the figure legend disagree with Table S3.

# Overlap match: a parsed CDS is assigned an annotation when the two intervals
# share at least half of the shorter one. Robust to the two gene callers
# disagreeing by a few codons at the ends.
apply_ann <- function(mapa, ann) {
  mapa$label <- ""
  mapa$cat   <- NA_character_
  hits <- integer(0)
  for (k in seq_len(nrow(ann))) {
    ov <- pmax(0, pmin(mapa$end, ann$end[k]) - pmax(mapa$start, ann$start[k]))
    len <- pmin(mapa$end - mapa$start, ann$end[k] - ann$start[k])
    i <- which(ov > 0.5 * len)
    if (length(i) == 0) {
      warning("no CDS overlaps ", ann$start[k], "-", ann$end[k],
              " (", ann$label[k], ")", call. = FALSE); next
    }
    i <- i[which.max(ov[i])]
    mapa$label[i] <- ann$label[k]; mapa$cat[i] <- ann$cat[k]
    hits <- c(hits, k)
  }
  attr(mapa, "matched") <- hits
  mapa
}

build <- function(path, ann) {
  m <- parse_gbk(path) %>% apply_ann(ann)
  m$Category <- relabel_map(ifelse(
    !is.na(m$cat), m$cat,
    ifelse(m$product != "hypothetical protein",
           "Structural / viral assembly", "Hypothetical / unannotated")))
  m$midpoint <- (m$start + m$end) / 2
  m
}
mapa_a <- build(path_a, ann_a) %>% asignar_niveles(min_dist = 3500)
mapa_b <- build(path_b, ann_b) %>% asignar_niveles(min_dist = 2500)

# ------------------------------------------------------------------------------
# 4. PANELS
# ------------------------------------------------------------------------------
base_panel <- function(m, brk, lims, title) {
  ggplot(m, aes(xmin = start, xmax = end, y = 1,
                fill = Category, forward = direccion_hebra)) +
    geom_gene_arrow(arrow_body_height = unit(5, "mm"),
                    arrowhead_height = unit(6, "mm"),
                    linewidth = 0.2, color = "black") +
    geom_segment(data = subset(m, label != ""),
                 aes(x = midpoint, xend = midpoint, y = 1.09, yend = label_y - 0.012),
                 inherit.aes = FALSE, colour = "grey55", linewidth = 0.25) +
    geom_text(data = subset(m, label != ""),
              aes(x = midpoint, y = label_y, label = label),
              inherit.aes = FALSE, hjust = 0.5, vjust = 0,
              size = SZ_GENE, fontface = "italic") +
    scale_x_continuous(breaks = brk, labels = function(x) paste0(x / 1000, " kb")) +
    coord_cartesian(xlim = lims, ylim = Y_MAP) +
    scale_fill_manual(values = paleta_global, drop = FALSE) +
    tema_isme + labs(x = NULL, title = title)
}

MT_L <- 17294; MT_R <- 21411          # the tandem methyltransferase array
pA <- base_panel(mapa_a, seq(0, 65000, 5000), c(0, 66500),
                 ttl6("A", NAME_a, "65.6 kb", "DE", "VIR", "culture-independent")) +
  # A statement of ABSENCE cannot be drawn as a bracket, so it is written out.
  annotate("segment", x = MT_L, xend = MT_R, y = 0.86, yend = 0.86,
           linewidth = 0.5, colour = "grey30") +
  annotate("segment", x = c(MT_L, MT_R), xend = c(MT_L, MT_R),
           y = 0.83, yend = 0.89, linewidth = 0.5, colour = "grey30") +
  annotate("text", x = (MT_L + MT_R) / 2, y = 0.79, vjust = 1,
           label = "three tandem methyltransferases \u00b7 no restriction endonuclease in 65.6 kb",
           size = SZ_NOTE, fontface = "italic", colour = "grey25")

GAJ_L <- 4947; GAJ_R <- 8860
pB <- base_panel(mapa_b, seq(0, 43000, 5000), c(0, 44000),
                 ttl6("B", NAME_b, "42.8 kb", "EL", "CUL", "culture-dependent", taxon = "Peduoviridae")) +
  # A bracket asserts that the system is complete, so it is drawn ONLY here:
  # L00454 is complete, confirmed, confidence high.
  annotate("segment", x = GAJ_L, xend = GAJ_R, y = 0.86, yend = 0.86, linewidth = 0.6) +
  annotate("segment", x = c(GAJ_L, GAJ_R), xend = c(GAJ_L, GAJ_R),
           y = 0.83, yend = 0.89, linewidth = 0.6) +
  annotate("text", x = (GAJ_L + GAJ_R) / 2, y = 0.79, vjust = 1,
           label = "Gabija \u00b7 complete (GajA + GajB) \u00b7 1 of 202 in the catalogue",
           size = SZ_NOTE, fontface = "italic", colour = "grey25")

# ------------------------------------------------------------------------------
# 5. VERIFICATION - read this before trusting the figure
# ------------------------------------------------------------------------------
check <- function(m, ann, lab) {
  cat("\n---", lab, "---\n")
  cat("CDS parsed from the GenBank file:", nrow(m), "\n")
  got <- attr(m, "matched")
  if (length(got) < nrow(ann))
    cat("*** UNMATCHED annotations:",
        paste(ann$start[-got], ann$end[-got], sep = "-", collapse = ", "), "***\n")
  print(as.data.frame(m %>% filter(!is.na(cat)) %>%
                        transmute(gene_num, start, end, strand,
                                  label = ifelse(label == "", "(unlabelled)", label),
                                  Category, product = substr(product, 1, 46))))
}
check(mapa_a, ann_a, "PANEL A - VIRFC_NODE_3")
check(mapa_b, ann_b, "PANEL B - EL-GLI_NODE_119")

# ------------------------------------------------------------------------------
# 6. LEGEND AND ASSEMBLY
# ------------------------------------------------------------------------------
p_legend <- ggplot(bind_rows(mapa_a, mapa_b),
                   aes(xmin = start, xmax = end, y = 1, fill = Category)) +
  geom_gene_arrow() + scale_fill_manual(values = paleta_global, drop = FALSE) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 15, face = "bold"),
        legend.key.size = unit(1.0, "cm"),
        legend.spacing.x = unit(0.8, "cm"),
        legend.box.margin = margin(t = 18, b = 18)) +
  guides(fill = guide_legend(ncol = 3, byrow = TRUE))

grab_legend <- function(p) {
  g <- try(cowplot::get_plot_component(p, "guide-box-bottom", return_all = TRUE),
           silent = TRUE)
  if (!inherits(g, "try-error") && !is.null(g) && !inherits(g, "zeroGrob")) return(g)
  cowplot::get_legend(p)
}
legend_final <- grab_legend(p_legend)
stopifnot(!is.null(legend_final))

figura_6 <- plot_grid(pA, pB, legend_final, ncol = 1, rel_heights = c(1, 1, 0.5))

# ------------------------------------------------------------------------------
# 7. EXPORT
# ------------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# EXPORT - the four formats ISME accepts, named as they are cited.
#
# ISME J, Author Guidelines > "Preparing your manuscript" > "Figures":
#   "Figures should be submitted in one of the following file formats:
#    .jpeg, .jpg, .png, .tiff, .pdf, or .eps."
#   "Save vector images as .eps or .pdf files and embed the fonts."
#   "Figure files should be named simply to match their citation
#    (e.g., fig1.tiff, fig2.eps)."
#
# Everything lands in ONE folder, ISME_figs/, inside the working directory,
# under the names the submission system expects: figN.pdf, .eps, .tiff, .png.
# Nothing above this block was changed, apart from removing setwd(MAPS_DIR) at the top.
#
# type = "cairo" on the raster devices. THIS IS THE FIX FOR THE BROKEN TIFF.
# ggsave() hands .tiff and .png to grDevices, which on Windows defaults to the
# GDI device; that device gives up on very large rasters and writes a blank or
# truncated file without erroring. Figure 2 is 16 x 14 in at 600 dpi = 9600 x
# 8400 px, which is where it breaks. cairo has no such limit, and it
# antialiases text the same way the PDF does.
#
# cairo_pdf, NOT pdf(): cairo EMBEDS the fonts, which is what the guideline
# asks for. The base device only references them, and it would also fail to
# render the bold / bold-italic mix in the bquote() titles.
#
# cairo_ps for EPS, NOT postscript(): the base one cannot embed fonts either.
#
# dpi is meaningless for PDF and EPS -- they have no pixels. Only width and
# height in inches fix the printed size, and both are unchanged.
# ---------------------------------------------------------------------------
OUT <- file.path(getwd(), "ISME_figs")
dir.create(OUT, showWarnings = FALSE)

export_isme <- function(plot, n, w, h) {
  if (!isTRUE(capabilities("cairo")))
    stop("This R build has no cairo support. Fonts cannot be embedded and the\n",
         "  large rasters will fail. Reinstall R with cairo before exporting.")
  b <- file.path(OUT, paste0("fig", n))
  ggsave(paste0(b, ".pdf"),  plot, width = w, height = h,
         device = cairo_pdf, bg = "white")
  ggsave(paste0(b, ".eps"),  plot, width = w, height = h,
         device = cairo_ps, fallback_resolution = 600, bg = "white")
  ggsave(paste0(b, ".tiff"), plot, width = w, height = h, dpi = 600,
         compression = "lzw", bg = "white", type = "cairo")
  ggsave(paste0(b, ".png"),  plot, width = w, height = h, dpi = 600,
         bg = "white", type = "cairo")
  cat("\n  ---------- fig", n, " ----------\n", sep = "")
  for (e in c(".pdf", ".eps", ".tiff", ".png"))
    cat(sprintf("    fig%d%-5s %8.2f MB\n", n, e,
                file.size(paste0(b, e)) / 1e6))
  cat("    folder: ", OUT, "\n", sep = "")
}

export_isme(figura_6, 6, 16, 8.6)

# Which palette_maps.R was actually loaded. Figures 4 and 5 read it from the
# WORKING DIRECTORY; Figure 6 reads its own copy from MAPS_DIR. They are not
# the same file, and that is why the three figures do not share a colour key.
# Check this line against the figure already in your manuscript BEFORE you
# upload anything.
pal_tag <- if ("Structural and other phage genes" %in% names(MAP_PALETTE))
             "NEW  (Structural and other phage genes / DGR retroelement)" else
             "OLD  (Structural / viral assembly / Hypermutation (DGR))"
cat("    palette loaded: ", pal_tag, "\n", sep = "")
