# =============================================================================
# 02_Fig_2.R — Figure 2.
#
# CHANGES (18/11)
#   1. Panel c: the PERMANOVA annotation now reads "Culture dependence"
#      instead of "Origin", matching the manuscript, where the term "origin"
#      has been removed throughout. The legend title is dropped altogether:
#      "Culture-dependent" and "Culture-independent" need no heading, and a
#      long title was being clipped at the left edge of the figure.
#   2. Host classes are harmonised to a single taxonomy before plotting
#      (section 0b). Host_Prediction merges three sources -- GTDB-Tk of the
#      source MAG, CHERRY-GTDB and CHERRY-NCBI -- and the three do not use
#      the same class names, so the same organism could appear under two
#      labels in panel H. 73 vOTUs whose order is Flavobacteriales were
#      labelled Bacteroidia (GTDB) while 910 identical ones were labelled
#      Flavobacteriia (NCBI).
#   3. New section 11: host-taxonomy audit, which re-checks the harmonised
#      column and should now report a single scheme.
#
# CHANGES (22/9)
#   1. Panel c: the PERMANOVA p-value is the exact one computed by 01_stats.R
#      (all 35 partitions; p = 0.029), and the "Origin | Zone" line is gone.
#   2. Sampling sites are labelled "Station" (legend of c, title of e1, x axis
#      of f and g), as in the text and the supplementary material.
#   3. Panel d prints "Odds ratio = 2.70" (spelled out, because OR is also the
#      code of the South Orkney Islands station), and its axis labels break
#      after the hyphen so they no longer touch.
#   4. Relative paths, and the input table is accepted with or without the
#      "_28_8" suffix: run from the folder that contains it.
#
# EARLIER CHANGES
#   1. No subtitles on any panel. Everything that was in a subtitle is in the
#      figure legend instead.
#   2. p_e1 draws all 15 intersections instead of the top 12. With 12 the
#      panel silently dropped EL+OR (58), EL+MB+OR (41) and EL+MB (8) --
#      107 vOTUs, 2.1% -- so the bars did not sum to n.
#
# CARRIED OVER FROM THE PREVIOUS TWO VERSIONS
#   - fmt_p prints the real p-value instead of flooring at 2.2e-16.
#   - p_d annotation on two lines so it fits inside the panel.
#   - p_e1 restricted to the 5,172 vOTUs assembled from CB22/GLI. Pooling all
#     11,100 gave 64.1% single-zone, but that is 81.6% for VIRFC against
#     47.0% for CB22, and the only reason is that VIRFC has no read set of
#     its own to be found in other zones with. Restricted: 47.8% / 24.5%.
#   - p_e2 is the original UpSet. Its set definitions are asymmetric
#     (BAC/CUL = read detection, VIR/LIS = assembly origin) and that goes in
#     the figure legend, not in the plot.
#   - This script computes no statistical test. Everything statistical is
#     read from stats_results.tsv.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(vegan); library(patchwork)
  library(ggupset);   library(stringr); library(RColorBrewer)
})
select <- dplyr::select

# Run from the folder that contains the input tables, e.g. setwd("path/to/R").
folder <- "Fig2"; dir.create(folder, showWarnings = FALSE)

# =============================================================================
# 0) READ DATA
# =============================================================================

# Input files. The scripts accept the names with or without the "_28_8" suffix,
# whichever exists in the working directory.
pick <- function(...) {
  f <- c(...); ok <- f[file.exists(f)]
  if (!length(ok)) stop("No file found: ", paste(f, collapse = " / "))
  cat("   input:", ok[1], "\n"); ok[1]
}
F_VOTU <- pick("07_votu_with_counts_28_8.tsv", "07_votu_with_counts.tsv")

votu_all <- read.delim(F_VOTU, stringsAsFactors = FALSE)
votu     <- votu_all %>% filter(Passes_catalog == "Yes")
stats    <- read.delim("stats_results.tsv", stringsAsFactors = FALSE)

# Only eight libraries were mapped against the catalogue: 4 CB22 and 4 GLI.
# VIRFC and LIS were never used as a read source, so cross-mapping onto their
# vOTUs is real signal (1,812 of 5,124 VIRFC and 627 of 804 LIS vOTUs are
# detected) but self-quantification does not exist for them.
SAMPLE_RE <- "^TPM_(DE|EL|MB|OR)_(CB22|GLI)$"
tpm_cols  <- grep(SAMPLE_RE, names(votu), value = TRUE)
stopifnot(length(tpm_cols) == 8)

votu <- votu %>%
  mutate(
    Life_Gen   = ifelse(grepl("^Temperate", Lifestyle), "Temperate", "Virulent"),
    Total_TPM  = rowSums(across(all_of(tpm_cols))),
    Viral_Fam  = str_extract(Taxonomy, "f_[^;]+"),
    Viral_Fam  = ifelse(is.na(Viral_Fam), "Unclassified", gsub("f_", "", Viral_Fam)),
    Host_Class = str_extract(Host_Prediction, "c__[^;]+"),
    Host_Class = ifelse(is.na(Host_Class) | grepl("Unknown", Host_Prediction),
                        "Unknown", gsub("c__", "", Host_Class))
  )

# =============================================================================
# 0b) HOST TAXONOMY HARMONISATION
#
#     Host_Prediction merges three sources: the GTDB-Tk taxonomy of the source
#     MAG, the GTDB-based CHERRY prediction and the NCBI-based one. The three
#     do not use the same class names, so without this step the same organism
#     can appear twice in panel H under two labels.
#
#     1,571 of the 2,157 assignments (72.8%) come from the NCBI branch, so NCBI
#     is the scheme with coverage and the GTDB- and MAG-derived minority is
#     translated into it. The case that matters is Flavobacteriales: GTDB puts
#     it inside Bacteroidia while NCBI gives it its own class, which was
#     splitting 73 vOTUs away from the 910 already labelled Flavobacteriia.
#
#     Classes whose circumscription differs between the two schemes are
#     resolved from the ORDER, which is the same in both. The rest are
#     one-to-one synonyms observed in this dataset.
# =============================================================================
ORDER_TO_NCBI <- c(
  "Flavobacteriales"   = "Flavobacteriia",
  "Bacteroidales"      = "Bacteroidia",
  "Cytophagales"       = "Cytophagia",
  "Chitinophagales"    = "Chitinophagia",
  "Sphingobacteriales" = "Sphingobacteriia",
  "Burkholderiales"    = "Betaproteobacteria")

GTDB_TO_NCBI <- c(
  "Actinomycetes"        = "Actinobacteria",
  "Cyanobacteriia"       = "Cyanophyceae",
  "Campylobacteria"      = "Epsilonproteobacteria",
  "Bacilli_A"            = "Mollicutes",
  "Desulfobacteria"      = "Deltaproteobacteria",
  "Desulfuromonadia"     = "Deltaproteobacteria",
  "Desulfotomaculia"     = "Clostridia",
  "Thermoanaerobacteria" = "Clostridia",
  "Leptospirae"          = "Spirochaetia",
  "Methanosarcinia"      = "Methanomicrobia")

host_class_raw <- votu$Host_Class

votu <- votu %>%
  mutate(
    Host_Order = str_extract(Host_Prediction, "o__[^;]+"),
    Host_Order = ifelse(is.na(Host_Order), NA_character_, sub("o__", "", Host_Order)),
    Host_Class = ifelse(
      Host_Source != "PhaBOX_NCBI" &
        Host_Class %in% c("Bacteroidia", "Gammaproteobacteria") &
        Host_Order %in% names(ORDER_TO_NCBI),
      unname(ORDER_TO_NCBI[Host_Order]), Host_Class),
    Host_Class = ifelse(
      Host_Source != "PhaBOX_NCBI" & Host_Class %in% names(GTDB_TO_NCBI),
      unname(GTDB_TO_NCBI[Host_Class]), Host_Class))

cat("\n  host classes relabelled: ", sum(host_class_raw != votu$Host_Class),
    " of ", sum(votu$Host_Class != "Unknown"), "\n", sep = "")
print(as.data.frame(tibble(before = host_class_raw, after = votu$Host_Class) %>%
                      filter(before != after) %>% count(before, after, sort = TRUE)),
      row.names = FALSE)

df_tpm <- votu %>%
  pivot_longer(all_of(tpm_cols), names_to = "Sample", values_to = "TPM") %>%
  filter(TPM > 0) %>%
  mutate(Zone_Sample  = sapply(strsplit(Sample, "_"), `[`, 2),
         Frac_Sample  = sapply(strsplit(Sample, "_"), `[`, 3),
         World_Mapped = ifelse(Frac_Sample == "CB22",
                               "Culture-independent", "Culture-dependent"),
         Sample_Clean = gsub("TPM_", "", Sample))

df_alpha <- df_tpm %>% group_by(Sample_Clean, World_Mapped) %>%
  summarise(Richness = n(), Shannon = diversity(TPM, index = "shannon"),
            .groups = "drop")

# --- print the p-value 01_stats.R computed, not R's 2.2e-16 printing floor -
SUP <- c("0"="⁰","1"="¹","2"="²","3"="³","4"="⁴",
         "5"="⁵","6"="⁶","7"="⁷","8"="⁸","9"="⁹",
         "-"="⁻")
sup <- function(s) paste0(SUP[strsplit(as.character(s), "")[[1]]], collapse = "")
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p >= 0.001) return(sprintf("p = %.3f", p))
  if (p == 0)   return("p < 1 × 10⁻³⁰⁰")
  e <- floor(log10(p)); m <- p / 10^e
  paste0("p = ", sprintf("%.1f", m), " × 10", sup(e))
}

# =============================================================================
# 1) THEME + COLOURS
# =============================================================================
BASE_SIZE <- 13
theme_pub <- theme_classic(base_size = BASE_SIZE) +
  theme(plot.title   = element_text(face = "bold", size = BASE_SIZE + 1,
                                    margin = margin(b = 6)),
        legend.title = element_text(size = BASE_SIZE - 1, face = "bold"),
        legend.text  = element_text(size = BASE_SIZE - 2),
        axis.title   = element_text(size = BASE_SIZE),
        axis.text    = element_text(size = BASE_SIZE - 1),
        strip.text   = element_text(size = BASE_SIZE - 1, face = "bold"))
source("palette_paper.R")

FRAC_LEVELS <- c("BAC", "CUL", "LIS", "VIR")

# =============================================================================
# 2) PANEL A
# =============================================================================
df_riq_frac <- votu %>% group_by(Fraction, Zone) %>%
  summarise(vOTUs = n(), .groups = "drop") %>%
  mutate(Fraction_lbl = factor(frac_lbl[Fraction], levels = FRAC_LEVELS))

p_a <- ggplot(df_riq_frac, aes(x = Fraction_lbl, y = vOTUs)) +
  geom_boxplot(aes(fill = Fraction), width = 0.4, alpha = 0.5,
               outlier.shape = NA, color = "grey50") +
  geom_jitter(aes(fill = Fraction), shape = 21, size = 3, stroke = 0.6,
              color = "black", width = 0.2) +
  scale_fill_manual(values = col_frac) +
  labs(title = "A | Assembled vOTUs", y = "vOTUs", x = "") +
  theme_pub + theme(legend.position = "none")

# =============================================================================
# 3) PANEL B
# =============================================================================
df_tpm_frac <- votu %>% group_by(Fraction, Zone) %>%
  summarise(TPM_sum = sum(Total_TPM), .groups = "drop") %>%
  mutate(Fraction_lbl = factor(frac_lbl[Fraction], levels = FRAC_LEVELS),
         logTPM = log10(TPM_sum + 1))

p_b <- ggplot(df_tpm_frac, aes(x = Fraction_lbl, y = logTPM)) +
  geom_boxplot(aes(fill = Fraction), width = 0.4, alpha = 0.5,
               outlier.shape = NA, color = "grey50") +
  geom_jitter(aes(fill = Fraction), shape = 21, size = 3, stroke = 0.6,
              color = "black", width = 0.2) +
  scale_fill_manual(values = col_frac) +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  labs(title = "B | Mapped abundance", y = "log10(CPM+1)", x = "") +
  theme_pub + theme(legend.position = "none")

# =============================================================================
# 4) PANEL C
#    "Origin" replaced by "Culture dependence" in the annotation, and the
#    legend title dropped: the two labels are self-explanatory and a long
#    title was clipped at the left edge.
# =============================================================================
otu_table <- df_tpm %>% select(Sample_Clean, vOTU_ID, TPM) %>%
  pivot_wider(names_from = vOTU_ID, values_from = TPM, values_fill = 0) %>%
  column_to_rownames("Sample_Clean")

pcoa <- cmdscale(vegdist(otu_table, method = "bray"), k = 2, eig = TRUE)
pcoa_coords <- as.data.frame(pcoa$points) %>% setNames(c("PCoA1", "PCoA2")) %>%
  rownames_to_column("Sample_Clean") %>%
  left_join(unique(df_tpm[, c("Sample_Clean", "World_Mapped", "Zone_Sample")]),
            by = "Sample_Clean") %>% rename(Zone = Zone_Sample)

perm_naive <- stats %>% filter(section == "beta diversity",
                               test == "PERMANOVA (bray)") %>% slice(1)

annot_label <- paste0(
  "PERMANOVA (Bray–Curtis)\n",
  "Culture dependence: R² = ", sprintf("%.3f", perm_naive$estimate), ", ",
  fmt_p(perm_naive$p))

p_c <- ggplot(pcoa_coords, aes(x = PCoA1, y = PCoA2, fill = World_Mapped)) +
  stat_ellipse(aes(color = World_Mapped), type = "t", linetype = 2,
               level = 0.95, show.legend = FALSE) +
  geom_point(aes(shape = Zone), size = 4, color = "black", stroke = 0.7) +
  scale_fill_manual(values = col_ori, name = NULL) +
  scale_color_manual(values = col_ori) +
  scale_shape_manual(values = c(21, 22, 23, 24), name = "Station") +
  annotate("text", x = -Inf, y = Inf, label = annot_label, hjust = -0.05,
           vjust = 1.1, size = BASE_SIZE * 0.28, fontface = "italic") +
  labs(title = "C | Community structure", x = "PCoA1", y = "PCoA2") +
  theme_pub +
  guides(fill = guide_legend(override.aes = list(shape = 21, color = "black"),
                             order = 1), shape = guide_legend(order = 2)) +
  theme(legend.position = "bottom", legend.box = "vertical")

# =============================================================================
# 5) PANEL D — annotation on two lines so it is not clipped
# =============================================================================
life_counts <- votu %>% count(Fraction_group, Life_Gen) %>%
  group_by(Fraction_group) %>% mutate(pct = n / sum(n)) %>% ungroup() %>%
  rename(World_Mapped = Fraction_group)

life_fisher <- stats %>% filter(section == "lifestyle",
                                comparison == "Culture-dependent vs Culture-independent") %>% slice(1)

p_d <- ggplot(life_counts, aes(x = World_Mapped, y = pct, fill = Life_Gen)) +
  geom_bar(stat = "identity", position = "fill", color = "black",
           linewidth = 0.3, width = 0.6) +
  scale_fill_manual(values = col_life, name = "") +
  scale_x_discrete(labels = function(x) sub("-", "-\n", x)) +
  scale_y_continuous(labels = scales::percent,
                     breaks = c(0, 0.25, 0.50, 0.75, 1.00)) +
  coord_cartesian(ylim = c(0, 1.24), clip = "off") +
  annotate("segment", x = 1, xend = 2, y = 1.04, yend = 1.04, linewidth = 0.5) +
  annotate("segment", x = 1, xend = 1, y = 1.00, yend = 1.04, linewidth = 0.5) +
  annotate("segment", x = 2, xend = 2, y = 1.00, yend = 1.04, linewidth = 0.5) +
  annotate("text", x = 1.5, y = 1.08,
           label = paste0("Odds ratio = ", sprintf("%.2f", life_fisher$estimate), "\n",
                          fmt_p(life_fisher$p)),
           size = BASE_SIZE * 0.27, fontface = "bold", vjust = 0,
           lineheight = 1.05) +
  labs(title = "D | Lifestyles", y = "vOTUs %", x = "") +
  theme_pub +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = BASE_SIZE - 3))

# =============================================================================
# 6) PANEL E1 — zone overlap. All 15 intersections (was 12).
#    Membership in zone Z = assembled at Z, or detected at Z by that
#    station's reads. Restricted to the 5,172 CB22/GLI-assembled vOTUs.
# =============================================================================
ZONES <- c("DE", "EL", "MB", "OR")
votu_cell <- votu %>% filter(Fraction %in% c("CB22", "GLI"))

zone_member <- function(df, z) {
  cols <- grep(paste0("^TPM_", z, "_"), tpm_cols, value = TRUE)
  df$Zone == z | rowSums(df[, cols, drop = FALSE]) > 0
}
zmat <- sapply(ZONES, function(z) zone_member(votu_cell, z))
colnames(zmat) <- ZONES

df_upset_zones <- tibble(
  vOTU_ID = votu_cell$vOTU_ID,
  Zones   = lapply(seq_len(nrow(zmat)), function(i) ZONES[zmat[i, ]]))
n_e1 <- nrow(df_upset_zones)

p_e1 <- ggplot(df_upset_zones, aes(x = Zones)) +
  geom_bar(fill = "grey75", color = "black", linewidth = 0.3) +
  scale_x_upset(n_intersections = 15) +
  labs(title = "E | Station overlap", y = "vOTUs", x = "") +
  theme_pub +
  theme_combmatrix(combmatrix.panel.point.color.fill = "black",
                   combmatrix.label.text = element_text(size = BASE_SIZE - 2))

# =============================================================================
# 7) PANEL E2 — fraction overlap. 11 intersections exist, all are drawn.
#    BAC = detected by CB22 reads; CUL = detected by GLI reads;
#    VIR = assembled from VIRFC; LIS = assembled from LIS.
#    VIR and LIS come from a single-valued column, so they are mutually
#    exclusive and a four-way intersection cannot exist. Do not report the
#    empty four-way intersection as a result.
# =============================================================================
frac_col <- function(f) grep(paste0("_", f, "$"), tpm_cols, value = TRUE)
set_BAC <- votu$vOTU_ID[rowSums(votu[, frac_col("CB22"), drop = FALSE]) > 0]
set_CUL <- votu$vOTU_ID[rowSums(votu[, frac_col("GLI"),  drop = FALSE]) > 0]
set_VIR <- votu$vOTU_ID[votu$Fraction == "VIRFC"]
set_LIS <- votu$vOTU_ID[votu$Fraction == "LIS"]

df_upset_fracs <- tibble(vOTU_ID = votu$vOTU_ID) %>% rowwise() %>%
  mutate(Fracs = list(as.character(na.omit(c(
    if (vOTU_ID %in% set_BAC) "BAC" else NA,
    if (vOTU_ID %in% set_CUL) "CUL" else NA,
    if (vOTU_ID %in% set_VIR) "VIR" else NA,
    if (vOTU_ID %in% set_LIS) "LIS" else NA))))) %>% ungroup()

p_e2 <- ggplot(df_upset_fracs, aes(x = Fracs)) +
  geom_bar(fill = "grey75", color = "black", linewidth = 0.3) +
  scale_x_upset(n_intersections = 12) +
  labs(title = "F | Fraction overlap", y = "vOTUs", x = "") +
  theme_pub +
  theme_combmatrix(combmatrix.panel.point.color.fill = "black",
                   combmatrix.label.text = element_text(size = BASE_SIZE - 2))

# =============================================================================
# 8) PANELS F/G
# =============================================================================
df_tax_vir_clean <- df_tpm %>% filter(Viral_Fam != "Unclassified")
top_vir <- names(sort(tapply(df_tax_vir_clean$TPM, df_tax_vir_clean$Viral_Fam,
                             sum), decreasing = TRUE)[1:10])
df_tax_vir_agg <- df_tax_vir_clean %>%
  mutate(Viral_Plot = factor(ifelse(Viral_Fam %in% top_vir, Viral_Fam, "Others"),
                             levels = c(top_vir, "Others"))) %>%
  group_by(Zone_Sample, World_Mapped, Life_Gen, Viral_Plot) %>%
  summarise(TPM_sum = sum(TPM), .groups = "drop") %>% rename(Zone = Zone_Sample)

p_f <- ggplot(df_tax_vir_agg, aes(x = Zone, y = TPM_sum, fill = Viral_Plot)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.3) +
  facet_grid(World_Mapped ~ Life_Gen) +
  scale_fill_manual(values = c(brewer.pal(10, "Paired"), "grey70"),
                    name = "Viral Family") +
  labs(title = "G | Viral taxonomy", y = "Relative abundance", x = "Station") +
  theme_pub + theme(axis.title.x = element_text(margin = margin(t = 6)),
                    panel.spacing = unit(0.5, "lines"))

df_tax_host_clean <- df_tpm %>% filter(Host_Class != "Unknown")
top_host <- names(sort(tapply(df_tax_host_clean$TPM, df_tax_host_clean$Host_Class,
                              sum), decreasing = TRUE)[1:10])
df_tax_host_agg <- df_tax_host_clean %>%
  mutate(Host_Plot = factor(ifelse(Host_Class %in% top_host, Host_Class, "Others"),
                            levels = c(top_host, "Others"))) %>%
  group_by(Zone_Sample, World_Mapped, Life_Gen, Host_Plot) %>%
  summarise(TPM_sum = sum(TPM), .groups = "drop") %>% rename(Zone = Zone_Sample)

p_g <- ggplot(df_tax_host_agg, aes(x = Zone, y = TPM_sum, fill = Host_Plot)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.3) +
  facet_grid(World_Mapped ~ Life_Gen) +
  scale_fill_manual(values = c(brewer.pal(10, "Set3"), "grey70"),
                    name = "Host Class") +
  labs(title = "H | Predicted host", y = "Relative abundance", x = "Station") +
  theme_pub + theme(axis.title.x = element_text(margin = margin(t = 6)),
                    panel.spacing = unit(0.5, "lines"))

# =============================================================================
# 9) ASSEMBLY + EXPORT
# =============================================================================
left_col  <- (p_a + p_b) / (p_c + p_d) / (p_e1 + p_e2) +
  plot_layout(heights = c(1, 1.2, 1))
right_col <- (p_f / p_g) + plot_layout(heights = c(1, 1)) &
  theme(plot.margin = margin(t = 5, b = 10, l = 5, r = 5))
final_plot <- (left_col | right_col) + plot_layout(widths = c(1.5, 1))

# ---------------------------------------------------------------------------
# EXPORT - PDF, EPS, SVG, TIFF, and PNG inside ISME_figs/.
# ---------------------------------------------------------------------------
OUT <- file.path(getwd(), "ISME_figs")
dir.create(OUT, showWarnings = FALSE)

export_isme <- function(plot, n, w, h) {
  if (!isTRUE(capabilities("cairo")))
    stop("This R build has no cairo support. Fonts cannot be embedded and the\n",
         "  large rasters will fail. Reinstall R with cairo before exporting.")
  b <- file.path(OUT, paste0("fig", n))
  svg_dev <- if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite else svg
  
  ggsave(paste0(b, ".pdf"),  plot, width = w, height = h,
         device = cairo_pdf, bg = "white")
  ggsave(paste0(b, ".eps"),  plot, width = w, height = h,
         device = cairo_ps, fallback_resolution = 600, bg = "white")
  ggsave(paste0(b, ".svg"),  plot, width = w, height = h,
         device = svg_dev, bg = "white")
  ggsave(paste0(b, ".tiff"), plot, width = w, height = h, dpi = 600,
         compression = "lzw", bg = "white", type = "cairo")
  ggsave(paste0(b, ".png"),  plot, width = w, height = h, dpi = 600,
         bg = "white", type = "cairo")
  
  cat("\n  ---------- fig", n, " ----------\n", sep = "")
  for (e in c(".pdf", ".eps", ".svg", ".tiff", ".png"))
    cat(sprintf("    fig%d%-5s %8.2f MB\n", n, e,
                file.size(paste0(b, e)) / 1e6))
  cat("    folder: ", OUT, "\n", sep = "")
}

export_isme(final_plot, 2, 16, 14)

# =============================================================================
# 10) DESCRIPTIVE REPORT — the numbers the text may quote for this figure
# =============================================================================
rule <- function(x) cat("\n", strrep("-", 74), "\n ", x, "\n",
                        strrep("-", 74), "\n", sep = "")

rule("1) CATALOGUE vs ABUNDANCE  (panels a, b)")
print(as.data.frame(votu %>% group_by(Fraction) %>%
                      summarise(vOTUs = n(), TPM_sum = sum(Total_TPM), .groups = "drop") %>%
                      mutate(vOTUs_pct = round(100 * vOTUs / sum(vOTUs), 2),
                             TPM_pct   = round(100 * TPM_sum / sum(TPM_sum), 2)) %>%
                      arrange(desc(vOTUs))))

rule("2) ALPHA DIVERSITY  (CB22 and GLI only, 4 vs 4)")
print(as.data.frame(df_alpha %>% group_by(World_Mapped) %>%
                      summarise(n = n(), Richness_mean = round(mean(Richness), 1),
                                Richness_sd = round(sd(Richness), 1),
                                Shannon_mean = round(mean(Shannon), 2),
                                Shannon_sd = round(sd(Shannon), 2), .groups = "drop")))

rule("3) LIFESTYLE  (panel d, all 11,100)")
print(as.data.frame(life_counts))
cat("Fisher: OR = ", round(life_fisher$estimate, 3), ", ",
    fmt_p(life_fisher$p), "\n", sep = "")

rule("4) BETA DIVERSITY  (read from stats_results.tsv)")
print(as.data.frame(stats %>% filter(section == "beta diversity")))

rule("5) PANEL E1")
e1_tab <- df_upset_zones %>% mutate(n_zones = lengths(Zones)) %>%
  count(n_zones) %>% mutate(pct = round(100 * n / sum(n), 2))
print(as.data.frame(e1_tab))
cat(e1_tab$pct[e1_tab$n_zones == 1], "% in one zone; ",
    e1_tab$pct[e1_tab$n_zones == 4], "% in all four; n = ", n_e1, "\n", sep = "")
cat("\nSingle-zone vOTUs, by station:\n")
print(as.data.frame(df_upset_zones %>% filter(lengths(Zones) == 1) %>%
                      mutate(Station = unlist(Zones)) %>% count(Station, sort = TRUE)))

rule("6) PANEL E2")
print(as.data.frame(df_upset_fracs %>% mutate(n_frac = lengths(Fracs)) %>%
                      count(n_frac) %>% mutate(pct = round(100 * n / sum(n), 2))))
det <- votu %>% mutate(
  in_cb = rowSums(across(all_of(frac_col("CB22")))) > 0,
  in_gl = rowSums(across(all_of(frac_col("GLI"))))  > 0)
cat("\nVIRFC vOTUs detected by cellular reads: ",
    sum(det$Fraction == "VIRFC" & (det$in_cb | det$in_gl)), " of ",
    sum(det$Fraction == "VIRFC"), "\n", sep = "")
cat("LIS   vOTUs detected by cellular reads: ",
    sum(det$Fraction == "LIS"   & (det$in_cb | det$in_gl)), " of ",
    sum(det$Fraction == "LIS"), "\n", sep = "")
cat("CB22/GLI-assembled vOTUs seen in BOTH read pools: ",
    sum(det$Fraction %in% c("CB22", "GLI") & det$in_cb & det$in_gl), " of ",
    sum(det$Fraction %in% c("CB22", "GLI")), "\n", sep = "")

rule("7) TOP TAXONOMY / HOST BY CULTURE DEPENDENCE  (panels g, h)")
print(as.data.frame(df_tpm %>% filter(Viral_Fam != "Unclassified") %>%
                      group_by(World_Mapped, Viral_Fam) %>% summarise(TPM = sum(TPM), .groups="drop") %>%
                      arrange(World_Mapped, desc(TPM)) %>% group_by(World_Mapped) %>% slice_head(n = 3) %>%
                      mutate(Pct = round(100 * TPM / sum(TPM), 2))))
print(as.data.frame(df_tpm %>% filter(Host_Class != "Unknown") %>%
                      group_by(World_Mapped, Host_Class) %>% summarise(TPM = sum(TPM), .groups="drop") %>%
                      arrange(World_Mapped, desc(TPM)) %>% group_by(World_Mapped) %>% slice_head(n = 3) %>%
                      mutate(Pct = round(100 * TPM / sum(TPM), 2))))

# =============================================================================
# 11) HOST TAXONOMY AUDIT  (panel H)
#
#     GTDB R226 and NCBI do not use the same class names. Betaproteobacteria,
#     Flavobacteriia, Cyanophyceae and Actinomycetes are NCBI classes that do
#     not exist in GTDB; Cyanobacteriia and Actinomycetia are their GTDB
#     counterparts, and GTDB places Flavobacteriales inside Bacteroidia.
#
#     If names from both schemes end up in the same column, one organism can
#     be counted twice under two different class labels and the proportions in
#     panel H are distorted. This block prints every class with its vOTU count
#     and says which scheme the names belong to.
# =============================================================================
rule("8) HOST TAXONOMY AUDIT  (panel H)")

cat("Host-related columns in the input table:\n  ")
cat(paste(grep("host", names(votu_all), ignore.case = TRUE, value = TRUE),
          collapse = ", "), "\n")

host_tab <- votu %>%
  filter(Host_Class != "Unknown") %>%
  count(Host_Class, sort = TRUE) %>%
  mutate(pct = round(100 * n / sum(n), 2))

cat("\nEvery host class, with vOTU counts (not just the ten plotted):\n")
print(as.data.frame(host_tab), row.names = FALSE)
cat("\n  distinct classes: ", nrow(host_tab),
    "  |  vOTUs with a class: ", sum(host_tab$n),
    "  |  without: ", sum(votu$Host_Class == "Unknown"), "\n", sep = "")

NCBI_ONLY <- c("Betaproteobacteria", "Deltaproteobacteria", "Epsilonproteobacteria",
               "Flavobacteriia", "Cyanophyceae", "Actinobacteria", "Mollicutes",
               "Sphingobacteriia", "Cytophagia", "Chitinophagia", "Saprospiria",
               "Methanomicrobia", "Erysipelotrichia", "Thermoprotei")
GTDB_ONLY <- c("Cyanobacteriia", "Actinomycetia", "Actinomycetes", "Campylobacteria",
               "Bacilli_A", "Leptospirae", "Desulfobacteria", "Desulfuromonadia",
               "Desulfotomaculia", "Thermoanaerobacteria", "Methanosarcinia",
               "Alphaproteobacteria_A", "Gammaproteobacteria_A", "Vampirovibrionia")

found_ncbi <- intersect(host_tab$Host_Class, NCBI_ONLY)
found_gtdb <- intersect(host_tab$Host_Class, GTDB_ONLY)

cat("\n  NCBI-only names found: ",
    if (length(found_ncbi)) paste(found_ncbi, collapse = ", ") else "none", "\n",
    "  GTDB-only names found: ",
    if (length(found_gtdb)) paste(found_gtdb, collapse = ", ") else "none", "\n",
    sep = "")

if (length(found_ncbi) && length(found_gtdb)) {
  cat("\n  >>> MIXED TAXONOMY STILL PRESENT after section 0b. Add the missing\n",
      "      names to ORDER_TO_NCBI or GTDB_TO_NCBI and re-run.\n", sep = "")
} else if (length(found_ncbi)) {
  cat("\n  >>> Host classes look NCBI-based. Say so in the Fig. 2 legend.\n", sep = "")
} else if (length(found_gtdb)) {
  cat("\n  >>> Host classes look GTDB-based. Say so in the Fig. 2 legend.\n", sep = "")
} else {
  cat("\n  >>> No diagnostic name present. Inspect Host_Prediction by hand.\n", sep = "")
}

if (all(c("Bacteroidia", "Flavobacteriia") %in% host_tab$Host_Class)) {
  cat("\n  WARNING: Bacteroidia and Flavobacteriia are both present. GTDB puts\n",
      "  Flavobacteriales inside Bacteroidia, so these two bars are very likely\n",
      "  the same organisms split across two taxonomies.\n", sep = "")
}

cat("\nFive raw Host_Prediction strings, to see the format:\n")
print(head(unique(votu$Host_Prediction[votu$Host_Class != "Unknown"]), 5))

cat("\nFig2: the five files are in ", OUT, "\n", sep = "")
