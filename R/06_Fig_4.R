# =============================================================================
# 07_Fig4.R — Figure 4. Genome maps of three virulent Caudoviricetes.
# Successor of 04_Fig_4_27_8.R. Ready to run.
#
# CHANGES 23/9
#   * orn and the two dcm genes take the category Data S2 gives them
#     ("Other / uncharacterised systems"); see the note next to the case_when.
#   * The panel titles give the vOTU length from Data S1 (126.7 kb), not the
#     length of the contig (127.2 kb).
#   * rfbA keeps the carbohydrate colour on purpose: rfbA, rfbE and wecB are
#     one sugar-biosynthesis operon, and the figure legend states that the
#     KEGG category of rfbA in Data S2 is biosynthesis of other secondary
#     metabolites. AMG colours follow KEGG; defence colours follow the tables.
#   * MAPS_DIR: the absolute paths live in one variable now.
#
# WHAT CHANGED EARLIER, AND WHY
#
#   1. THE CONTEXT BARS NOW SAY WHICH GENOME THEY BELONG TO.
#      The bars were already built correctly -- plot_grid(overview_B, pB) puts
#      B's overview above B's panel -- but they are drawn in the gap between
#      the previous panel's axis and their own title, so the eye assigns them
#      upwards. A reader saw "Region shown in b" under panel a and concluded
#      that a contains that region. It does not: a is a 29.2 kb contig and
#      that bar represents a 101.7 kb one. Naming the contig on the bar makes
#      the misreading impossible regardless of spacing.
#
#   2. ASYMMETRIC MARGINS on the overview bars, so each groups visually with
#      the panel below it instead of floating between two panels.
#
#   3. TITLES: contig name first, one set of parentheses. Panel a had two.
#
#   4. get_legend() rewritten. The old one matched a grob called "guide-box",
#      which ggplot2 3.5 renamed; it returns nothing on current versions and
#      the figure loses its legend without any error.
#
#   Colours and every gene assignment are untouched. Category names come from
#   palette_maps.R, shared with Figures 5 and 6.
#
# THE GENE ASSIGNMENTS ARE MANUAL AND THAT IS DELIBERATE
#   Categories are given per gene_num because GenBank product names alone do
#   not resolve them. The folE/queE/queD group in panel a is a case in point:
#   DefenseFinder calls those three a Dpd system, and the core-function veto
#   reclassified them as AMGs because they are the shared preQ0 biosynthesis
#   pathway (Table S5, "Corrected during audit"). The figure shows the result
#   of that curation, so it must not be re-derived from the raw annotation.
# =============================================================================

library(ggplot2)
library(dplyr)
library(gggenes)
library(cowplot)

# ------------------------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------------------------
# One line to change for the repository: MAPS_DIR is the only absolute path.
MAPS_DIR <- "C:/Users/gasch/Documents/R/vOTUs-R/MAPAS"

path_A <- file.path(MAPS_DIR, "liticos/MCHA-DE-14-CB22_NODE_1150_length_29199_cov_7.788945/vOTU_viral_15591.gbf")
path_B <- file.path(MAPS_DIR, "liticos/MCHA-DE-14-VIRFC_NODE_1_length_101954_cov_30.149727/vOTU_viral_10366.gbf")
path_C <- file.path(MAPS_DIR, "liticos/MCHA-MB-1-CB22_NODE_23_length_127234_cov_9.834242/vOTU_viral_14385.gbf")

NAME_A <- "MCHA-DE-14-CB22_NODE_1150"
NAME_B <- "MCHA-DE-14-VIRFC_NODE_1"
NAME_C <- "MCHA-MB-1-CB22_NODE_23"

# ------------------------------------------------------------------------------
# 2. PALETTE AND PARSER
# ------------------------------------------------------------------------------
source("palette_maps.R")
paleta_global <- MAP_PALETTE

parse_gbk <- function(file) {
  lines <- readLines(file)
  cds_idx <- grep("^     CDS ", lines)
  data <- do.call(rbind, lapply(cds_idx, function(i) {
    loc_line <- lines[i]
    nums <- as.numeric(regmatches(loc_line, gregexpr("[0-9]+", loc_line))[[1]])
    prod <- "hypothetical protein"
    for (j in (i + 1):(i + 20)) {
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

# Figure 5 is 16 x 9.5 in and this one was 16 x 16, so at a common print
# width every glyph here came out a third smaller than the same glyph there.
# The height is now 13 and the type is bigger; see also the tighter ylim in
# the map panels, which was wasting the top of every one of them.
SZ_GENE  <- 5.0    # italic gene labels on the arrows
SZ_AXIS  <- 13     # kb tick labels
SZ_TITLE <- 16     # panel titles
SZ_OVER_NAME <- 4.4   # contig name over the context bar
SZ_OVER_LAB  <- 4.1   # region captions under it
Y_MAP    <- c(0.90, 1.45)   # was c(0.85, 1.60): nothing is drawn above 1.35

tema_isme <- theme_minimal() + theme(
  panel.grid = element_blank(), axis.line.x = element_line(color = "black"),
  axis.text.y = element_blank(), axis.title = element_blank(),
  # The kb tick labels were rendering at the theme default and came out too
  # small to read next to 4.5-pt gene labels.
  axis.text.x = element_text(size = SZ_AXIS, colour = "grey20"),
  axis.ticks.x = element_line(colour = "grey20", linewidth = 0.4),
  legend.position = "none",
  plot.title = element_text(face = "bold", size = SZ_TITLE),
  plot.margin = margin(t = 10, r = 20, b = 20, l = 20)
)

# The overview bar belongs to the panel BELOW it. Space above, none below.
tema_overview <- theme_void() +
  theme(plot.margin = margin(t = 26, r = 20, b = 0, l = 20))

# The context bar went through two bad versions. The first was a 0.10-tall
# rectangle in a 0.65 y range inside a panel given 22% of its row: a hairline
# on paper. The second was thicker but put the contig name and the region
# captions both ABOVE the bar, where they collided with each other and with
# the highlighted block. Here the name sits above and the region captions
# below, so nothing can overlap whatever the genome length or region width.
# No alpha anywhere. cairo_ps has to rasterise any semi-transparent object and
# falls back to a bitmap at fallback_resolution for that region, so a single
# alpha = 0.85 rectangle turns part of a vector EPS into pixels. #6D6D6D is
# what #595959 at 85% over grey88 resolves to, so the figure looks identical
# and the EPS stays fully vector.
overview_bar <- function(genome_len, name, kb_label, regions, labels) {
  p <- ggplot() +
    annotate("rect", xmin = 0, xmax = genome_len, ymin = 1.10, ymax = 1.32,
             fill = "grey88", color = "black", linewidth = 0.5)
  for (i in seq_along(regions)) {
    p <- p + annotate("rect", xmin = regions[[i]][1], xmax = regions[[i]][2],
                      ymin = 1.10, ymax = 1.32, fill = "#6D6D6D",
                      color = "black", linewidth = 0.6)
  }
  for (i in seq_along(labels)) {
    p <- p + annotate("text", x = labels[[i]]$x, y = 0.88, label = labels[[i]]$text,
                      size = SZ_OVER_LAB, fontface = "italic", vjust = 1)
  }
  p +
    annotate("text", x = 0, y = 1.52, hjust = 0, vjust = 0,
             label = paste0(name, "  \u00b7  ", kb_label),
             size = SZ_OVER_NAME, fontface = "bold") +
    # theme_void() has no axes, so anything outside ylim is silently clipped.
    # After the type got bigger and the figure shorter, the contig name and the
    # region captions were being cut off at both ends. The range is now wide
    # enough that neither can reach an edge.
    coord_cartesian(xlim = c(0, genome_len), ylim = c(0.35, 2.10)) +
    tema_overview
}

# THE FRACTION IS ADDED TO THE TITLE, NOT SUBSTITUTED FOR THE ORIGIN.
# "Cult.-indep." was already correct on all four panels: BAC (CB22) and VIR
# (VIRFC) are both culture-independent, and all three contigs here come from
# one or the other. What the four identical labels did not say is WHICH
# fraction, and the fraction is the variable Figure 3e is built on -- BAC sits
# at 29.2% of loci on the exclusion axis and VIR at 10.3%. Both labels are
# therefore kept: the origin, which is what the lifestyle analysis uses, and
# the fraction, which is what the axis analysis uses. They are different
# variables and neither replaces the other.
ttl <- function(letter, name, kb, zone, frac, extra = NULL) {
  tail <- paste0(" (", kb, " \u00b7 ", zone, " \u00b7 ", frac,
                 " \u00b7 culture-independent)")
  if (!is.null(extra)) tail <- paste0(tail, "  |  ", extra)
  bquote(bold(.(paste0(letter, " | Virulent "))) *
           bolditalic("Caudoviricetes") *
           bold(.(paste0(" \u00b7 ", name, tail))))
}

# ------------------------------------------------------------------------------
# 3. MAPS
# ------------------------------------------------------------------------------

# --- vOTU A (29.2 kb) ---------------------------------------------------------
mapa_A <- parse_gbk(path_A) %>% mutate(
  Category = relabel_map(case_when(
    gene_num %in% c(23, 24, 25)       ~ "AMG: Metabolism of cofactors and vitamins",
    product != "hypothetical protein" ~ "Structural / viral assembly",
    TRUE                              ~ "Hypothetical / unannotated")),
  label = case_when(gene_num == 23 ~ "folE", gene_num == 24 ~ "queE",
                    gene_num == 25 ~ "queD", TRUE ~ ""),
  y_pos = 1.15
)
pA <- ggplot(mapa_A, aes(xmin = start, xmax = end, y = 1,
                         fill = Category, forward = direccion_hebra)) +
  geom_gene_arrow(arrow_body_height = unit(5, "mm"),
                  arrowhead_height = unit(6, "mm"), linewidth = 0.2, color = "black") +
  geom_text(aes(x = (start + end) / 2, y = y_pos, label = label),
            angle = 45, hjust = 0, size = SZ_GENE, fontface = "italic") +
  scale_x_continuous(breaks = seq(0, 30000, 5000),
                     labels = function(x) paste0(x / 1000, " kb")) +
  coord_cartesian(xlim = c(0, 32000), ylim = Y_MAP) +
  scale_fill_manual(values = paleta_global, drop = FALSE) +
  tema_isme + labs(title = ttl("a", NAME_A, "29.2 kb", "DE", "BAC"))

# --- vOTU B (101.7 kb) --------------------------------------------------------
mapa_B <- parse_gbk(path_B) %>% mutate(
  Category = relabel_map(case_when(
    gene_num == 56                    ~ "AMG: Metabolism of cofactors and vitamins",
    gene_num %in% c(51, 57, 60)       ~ "AMG: Carbohydrate metabolism",
    product != "hypothetical protein" ~ "Structural / viral assembly",
    TRUE                              ~ "Hypothetical / unannotated")),
  label = case_when(gene_num == 51 ~ "rfbA", gene_num == 56 ~ "ubiG",
                    gene_num == 57 ~ "rfbE", gene_num == 60 ~ "wecB", TRUE ~ ""),
  y_pos = 1.15
)
zoom_B_start <- min(mapa_B$start[mapa_B$gene_num == 48]) - 500
zoom_B_end   <- max(mapa_B$end[mapa_B$gene_num == 62]) + 500

pB <- ggplot(mapa_B, aes(xmin = start, xmax = end, y = 1,
                         fill = Category, forward = direccion_hebra)) +
  geom_gene_arrow(arrow_body_height = unit(5, "mm"),
                  arrowhead_height = unit(6, "mm"), linewidth = 0.2, color = "black") +
  geom_text(aes(x = (start + end) / 2, y = y_pos, label = label),
            angle = 45, hjust = 0, size = SZ_GENE, fontface = "italic") +
  scale_x_continuous(breaks = seq(zoom_B_start, zoom_B_end, 2000),
                     labels = function(x) paste0(round(x / 1000, 1), " kb")) +
  coord_cartesian(xlim = c(zoom_B_start, zoom_B_end), ylim = Y_MAP) +
  scale_fill_manual(values = paleta_global, drop = FALSE) +
  tema_isme + labs(title = ttl("b", NAME_B, "101.7 kb", "DE", "VIR", "zoom"))

overview_B <- overview_bar(
  genome_len = max(mapa_B$end), name = NAME_B, kb_label = "101.7 kb",
  regions = list(c(zoom_B_start, zoom_B_end)),
  labels  = list(list(x = mean(c(zoom_B_start, zoom_B_end)),
                      text = "region shown in b, below")))
pB_full <- plot_grid(overview_B, pB, ncol = 1, rel_heights = c(0.46, 1))

# --- vOTU C (126.7 kb) --------------------------------------------------------
mapa_C <- parse_gbk(path_C) %>% mutate(
  Category = relabel_map(case_when(
    # CATEGORIES FOLLOW Data S2 (23/9). orn is a DEDDh exonuclease called by
    # DefenseFinder with no system and no locus, so Data S2 files it under
    # "Other / uncharacterised systems", not nucleic acid degradation. The dcm
    # at 35.8 kb and the one at 94.2 kb are Druantia_II calls whose loci
    # (L00878, L00881) are 1-of-2 fragments, so they keep their gene-level
    # category as well. Genes 53, 105 and 113 stay in genome modification
    # because they sit inside COMPLETE orphan-methyltransferase loci (L00877,
    # L00879/L00880), which is the exception declared in Supplementary
    # Methods S15.
    gene_num == 3                          ~ "Defence: Other / uncharacterised",
    gene_num %in% c(84, 145)               ~ "Defence: Other / uncharacterised",
    gene_num %in% c(53, 105, 113)          ~ "Defense: Genome modification",
    gene_num == 112                        ~ "AMG: Metabolism of cofactors and vitamins",
    product != "hypothetical protein"      ~ "Structural / viral assembly",
    TRUE                                   ~ "Hypothetical / unannotated")),
  label = case_when(
    gene_num == 3   ~ "orn",      gene_num == 53  ~ "N6-MTase",
    gene_num == 84  ~ "dcm",      gene_num == 105 ~ "MTase",
    gene_num == 112 ~ "cobC",     gene_num == 113 ~ "N4-MTase",
    gene_num == 145 ~ "dcm",      TRUE            ~ ""),
  y_pos = 1.15,
  x_text = case_when(
    gene_num == 112 ~ ((start + end) / 2) - 900,   # cobC pushed left
    gene_num == 113 ~ ((start + end) / 2) + 900,   # N4-MTase pushed right
    TRUE            ~ (start + end) / 2)
)

map_C_panel <- function(brk, lims, letter, locus) {
  ggplot(mapa_C, aes(xmin = start, xmax = end, y = 1,
                     fill = Category, forward = direccion_hebra)) +
    geom_gene_arrow(arrow_body_height = unit(5, "mm"),
                    arrowhead_height = unit(6, "mm"), linewidth = 0.2, color = "black") +
    geom_text(aes(x = x_text, y = y_pos, label = label),
              angle = 45, hjust = 0, size = SZ_GENE, fontface = "italic") +
    scale_x_continuous(breaks = brk, labels = function(x) paste0(x / 1000, " kb")) +
    coord_cartesian(xlim = lims, ylim = Y_MAP) +
    scale_fill_manual(values = paleta_global, drop = FALSE) +
    tema_isme + labs(title = ttl(letter, NAME_C, "126.7 kb", "MB", "BAC", locus))
}
pC1 <- map_C_panel(seq(0, 55000, 10000),     c(0, 56000),      "c", "locus 1")
pC2 <- map_C_panel(seq(90000, 100000, 2000), c(90000, 100000), "d", "locus 2")

genome_C_len <- max(mapa_C$end)
overview_C <- overview_bar(
  genome_len = genome_C_len, name = NAME_C, kb_label = "126.7 kb",
  regions = list(c(0, 56000), c(90000, 100000)),
  labels  = list(list(x = 28000, text = "locus 1, in c below"),
                 list(x = 95000, text = "locus 2, in d below")))
# Panel d carries a single label against panel c's seven, so giving them equal
# height left d as mostly blank paper.
pC_full <- plot_grid(overview_C, pC1, pC2, ncol = 1, rel_heights = c(0.44, 1.15, 0.75))

# ------------------------------------------------------------------------------
# 4. LEGEND
# ------------------------------------------------------------------------------
all_maps_combined <- bind_rows(mapa_A, mapa_B, mapa_C)
p_legend <- ggplot(all_maps_combined, aes(xmin = start, xmax = end, y = 1, fill = Category)) +
  geom_gene_arrow() + scale_fill_manual(values = paleta_global, drop = FALSE) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        # 12 pt read fine on a 9.5-inch figure and small on this 16-inch one.
        legend.text = element_text(size = 15, face = "bold"),
        legend.key.size = unit(1.0, "cm"),
        legend.spacing.x = unit(1.0, "cm"), legend.spacing.y = unit(0.5, "cm"),
        legend.box.margin = margin(t = 25, b = 25, l = 10, r = 10)) +
  guides(fill = guide_legend(ncol = 3, byrow = TRUE))

# ggplot2 3.5 renamed the legend grob from "guide-box" to "guide-box-bottom",
# so the old hand-rolled extractor silently returned nothing and the figure
# came out with no legend and no error. Try the current accessor first.
grab_legend <- function(p) {
  g <- try(cowplot::get_plot_component(p, "guide-box-bottom", return_all = TRUE),
           silent = TRUE)
  if (!inherits(g, "try-error") && !is.null(g) && !inherits(g, "zeroGrob")) return(g)
  cowplot::get_legend(p)
}
legend_final <- grab_legend(p_legend)
stopifnot(!is.null(legend_final))

# ------------------------------------------------------------------------------
# 5. ASSEMBLE AND SAVE
# ------------------------------------------------------------------------------
figura_4_mapas <- plot_grid(pA, pB_full, pC_full, legend_final,
                            ncol = 1, rel_heights = c(1, 1.55, 2.35, 1.05))

# ------------------------------------------------------------------------------
# EXPORT — PNG for looking at, TIFF and EPS for the journal
#
# EPS via cairo_ps rather than the base postscript device, because the base one
# cannot embed fonts or handle any transparency and silently drops what it
# cannot draw. cairo_ps is a build-time option in R, so it is checked rather
# than assumed: on a build without it the EPS is skipped with a message instead
# of erroring after the other two files are already written.
#
# fallback_resolution only matters for objects cairo cannot express as vectors
# (semi-transparent fills). There are none left in this figure, so nothing here
# should rasterise; it is set anyway so that adding an alpha later degrades to
# a 600-dpi bitmap rather than to a blank area.
# ------------------------------------------------------------------------------
save_all <- function(plot, base, w, h) {
  ggsave(paste0(base, ".png"),  plot, width = w, height = h, dpi = 600, bg = "white")
  ggsave(paste0(base, ".tiff"), plot, width = w, height = h, dpi = 600,
         compression = "lzw", bg = "white")
  if (isTRUE(capabilities("cairo"))) {
    ggsave(paste0(base, ".eps"), plot, width = w, height = h,
           device = cairo_ps, fallback_resolution = 600, bg = "white")
    cat("  wrote ", base, ".png / .tiff / .eps\n", sep = "")
  } else {
    warning("This R build has no cairo support, so no EPS was written. ",
            "Install R with cairo, or open the TIFF in Illustrator/Inkscape ",
            "and export EPS from there.", call. = FALSE)
    cat("  wrote ", base, ".png / .tiff  (EPS skipped: no cairo)\n", sep = "")
  }
}

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
save_all(figura_4_mapas, paste0("Figure_4_GenomicMaps_", stamp), 16, 13)
save_all(figura_4_mapas, "Figure_4_GenomicMaps", 16, 13)
cat("Send THIS file: Figure_4_GenomicMaps_", stamp, ".png\n", sep = "")

# ------------------------------------------------------------------------------
# 6. WHAT IS ON EACH PANEL — printed so the legend text can be checked
# ------------------------------------------------------------------------------
for (nm in c("mapa_A", "mapa_B", "mapa_C")) {
  cat("\n--", nm, "--\n")
  print(as.data.frame(get(nm) %>% count(Category, .drop = FALSE) %>% filter(n > 0)))
  lab <- get(nm) %>% filter(label != "") %>% select(gene_num, label, Category)
  if (nrow(lab)) print(as.data.frame(lab))
}
