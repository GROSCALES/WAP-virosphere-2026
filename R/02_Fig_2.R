# =============================================================================
# 02_Fig_2.R — Figure 2.
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
SUP <- c("0"="\u2070","1"="\u00b9","2"="\u00b2","3"="\u00b3","4"="\u2074",
         "5"="\u2075","6"="\u2076","7"="\u2077","8"="\u2078","9"="\u2079",
         "-"="\u207b")
sup <- function(s) paste0(SUP[strsplit(as.character(s), "")[[1]]], collapse = "")
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p >= 0.001) return(sprintf("p = %.3f", p))
  if (p == 0)   return("p < 1 \u00d7 10\u207b\u00b3\u2070\u2070")
  e <- floor(log10(p)); m <- p / 10^e
  paste0("p = ", sprintf("%.1f", m), " \u00d7 10", sup(e))
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
  "PERMANOVA (Bray\u2013Curtis)\n",
  "Origin: R\u00b2 = ", sprintf("%.3f", perm_naive$estimate), ", ", fmt_p(perm_naive$p))

p_c <- ggplot(pcoa_coords, aes(x = PCoA1, y = PCoA2, fill = World_Mapped)) +
  stat_ellipse(aes(color = World_Mapped), type = "t", linetype = 2,
               level = 0.95, show.legend = FALSE) +
  geom_point(aes(shape = Zone), size = 4, color = "black", stroke = 0.7) +
  scale_fill_manual(values = col_ori, name = "Origin") +
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
# Nothing above this block was changed.
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

# WARNING, FIGURE 2 ONLY. Panels A and B use alpha = 0.5 in geom_boxplot().
# EPS cannot carry transparency at all, so cairo_ps rasterises those two
# boxplot areas at 600 dpi; the rest of the EPS stays vector. The PDF and the
# TIFF are unaffected. Removing the alpha would make the EPS fully vector but
# would change how the figure looks, so it is NOT done here. Submit the PDF or
# the TIFF for this figure and the question does not arise.

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

rule("7) TOP TAXONOMY / HOST BY ORIGIN  (panels f, g)")
print(as.data.frame(df_tpm %>% filter(Viral_Fam != "Unclassified") %>%
                      group_by(World_Mapped, Viral_Fam) %>% summarise(TPM = sum(TPM), .groups="drop") %>%
                      arrange(World_Mapped, desc(TPM)) %>% group_by(World_Mapped) %>% slice_head(n = 3) %>%
                      mutate(Pct = round(100 * TPM / sum(TPM), 2))))
print(as.data.frame(df_tpm %>% filter(Host_Class != "Unknown") %>%
                      group_by(World_Mapped, Host_Class) %>% summarise(TPM = sum(TPM), .groups="drop") %>%
                      arrange(World_Mapped, desc(TPM)) %>% group_by(World_Mapped) %>% slice_head(n = 3) %>%
                      mutate(Pct = round(100 * TPM / sum(TPM), 2))))

cat("\nFig2: the four files are in ", OUT, "\n", sep = "")
