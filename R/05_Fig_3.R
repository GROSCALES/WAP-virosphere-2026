# =============================================================================
# 05_Fig_3.R — Figure 3. THE canonical script. Overwrite it, never branch it.
#
# RUN ORDER
#   01_stats.R  ->  03_stats_loci.R  ->  04_stats_fig3.R  ->  05_Fig_3.R
#
# FIVE PANELS
#   a  Density per Mb by category, mirrored by origin, split by lifestyle.
#      Defence counted as architecturally complete loci, AMGs as genes.
#   b  Odds ratios on the same rows, two columns: origin, and lifestyle.
#   c  The same contrast collapsed over categories.
#   d  System families of the 773 complete loci, faceted by axis.
#   e  Axis composition across fractions.
#
# VERSION 23/9. What changed over the 31/8 script:
#   1. Panel b now states that its odds ratios are STRATIFIED. Since
#      03_stats_loci.R and 04_stats_fig3.R switched to Mantel-Haenszel, the
#      column headers read "within lifestyles" and "within fractions" and the
#      x axis says Mantel-Haenszel, exactly as in Table S4a. A guard stops the
#      script if it is handed the old crude-Fisher files.
#   2. The two in-panel annotations are spelled out: "origin, given lifestyle"
#      in c and "lifestyle, given fraction" in e. The bar "|" read as a pipe.
#   3. No absolute path. Set the working directory to the folder with the
#      tables and source the file; input names are tried with and without the
#      "_28_8" suffix.
#
# ABUNDANCE IS NOT IN THIS FIGURE, AND THAT IS DELIBERATE
#   Two attempts were made and both dropped. A category x station heatmap
#   separated nothing that panel a does not already separate and cost a
#   quarter of the page. A CPM row under panel c pointed in exactly the same
#   direction as the per-Mb row for both classes, so it doubled that panel's
#   height to repeat a conclusion. Abundance lives in Table S4c.
#
# NAMING
#   The gene CLASS is Defence vs AMG. The two AXES are Self-protection vs
#   Exclusion. See palette_paper.R.
#
# CRISPR-Cas AND HOST MANIPULATION are absent from a, b and d because they
# have ZERO architecturally complete loci on the two antiviral axes. They
# exist as genes and as loci in Data S2 and Data S3. That absence is a
# result, stated in the legend.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(cowplot); library(scales)
})
select <- dplyr::select
source("palette_paper.R")
folder <- "Fig3"; dir.create(folder, showWarnings = FALSE)

need <- c("votu_n_defence_loci.tsv", "07_long_table_loci_for_fig3.tsv",
          "stats_results_loci_defence.tsv", "stats_results_fig3.tsv",
          "stats_results.tsv")
miss <- need[!file.exists(need)]
if (length(miss)) stop("Missing: ", paste(miss, collapse = ", "),
                       "\n  Run 01_stats.R, 03_stats_loci.R and 04_stats_fig3.R first.")

pick <- function(...) {
  f <- c(...); ok <- f[file.exists(f)]
  if (!length(ok)) stop("No file found: ", paste(f, collapse = " / "))
  cat("   input:", ok[1], "\n"); ok[1]
}
F_VOTU <- pick("07_votu_with_counts_28_8.tsv", "07_votu_with_counts.tsv")
F_LONG <- pick("07_long_table_FINAL_28_8.tsv", "07_long_table_FINAL.tsv")
F_LOCI <- pick("10_loci_defensa_28_8.tsv", "10_loci_defensa.tsv")

BASE <- 12
SZ_NUM   <- BASE * 0.27    # numbers printed on or beside bars
SZ_SMALL <- BASE * 0.22    # n = / OR = annotations
SZ_CATS  <- BASE - 1       # category names on the y axis

theme_pub <- theme_classic(base_size = BASE) +
  theme(legend.title = element_text(size = BASE - 2, face = "bold"),
        legend.text  = element_text(size = BASE - 3),
        axis.title   = element_text(size = BASE - 1),
        axis.text    = element_text(size = BASE - 1),
        strip.text   = element_text(size = BASE - 1, face = "bold"),
        strip.background = element_rect(fill = "grey93", colour = NA),
        plot.margin  = margin(t = 20, r = 6, b = 4, l = 4))

grab_legend <- function(p) {
  g <- try(cowplot::get_plot_component(p, "guide-box-bottom", return_all = TRUE),
           silent = TRUE)
  if (!inherits(g, "try-error") && !is.null(g) && !inherits(g, "zeroGrob")) return(g)
  cowplot::get_legend(p)
}
star <- function(q) ifelse(is.na(q), "", ifelse(q < 1e-3, "***",
                                                ifelse(q < 1e-2, "**", ifelse(q < 0.05, "*", ""))))

# =============================================================================
# READ
# =============================================================================
votu <- read.delim(F_VOTU, stringsAsFactors = FALSE) %>%
  filter(Passes_catalog == "Yes") %>%
  left_join(read.delim("votu_n_defence_loci.tsv", stringsAsFactors = FALSE),
            by = "vOTU_ID") %>%
  mutate(n_Defence_loci = replace_na(n_Defence_loci, 0),
         Origin = factor(ifelse(Fraction_group == "Culture-dependent",
                                "Culture-dependent", "Culture-independent"),
                         levels = c("Culture-independent", "Culture-dependent")),
         Life   = factor(ifelse(grepl("^Temperate", Lifestyle), "Temperate", "Virulent"),
                         levels = c("Temperate", "Virulent")),
         Frac   = factor(frac_lbl[Fraction], levels = frac_lev_cell))
tpm_cols <- grep("^TPM_(DE|EL|MB|OR)_(CB22|GLI)$", names(votu), value = TRUE)
stopifnot(length(tpm_cols) == 8)

long <- bind_rows(
  read.delim(F_LONG, stringsAsFactors = FALSE) %>%
    filter(Passes_catalog == "Yes", gene_class == "AVG") %>%
    select(vOTU_ID, category_final) %>% mutate(Class = "AMG"),
  read.delim("07_long_table_loci_for_fig3.tsv", stringsAsFactors = FALSE) %>%
    select(vOTU_ID, category_final) %>% mutate(Class = "Defence")) %>%
  left_join(votu %>% select(vOTU_ID, Origin, Life, Frac), by = "vOTU_ID") %>%
  mutate(Class = factor(Class, levels = c("Defence", "AMG")))

loci   <- read.delim(F_LOCI, stringsAsFactors = FALSE)
S_LOCI <- read.delim("stats_results_loci_defence.tsv", stringsAsFactors = FALSE)
S_FIG3 <- read.delim("stats_results_fig3.tsv", stringsAsFactors = FALSE)

# GUARD: the figure states Mantel-Haenszel in its legend, so it refuses to
# draw crude odds ratios. If this stops, re-run 03_stats_loci.R (23/9).
if (!any(grepl("Mantel", S_LOCI$test[S_LOCI$section == "category enrichment (loci)"])))
  stop("stats_results_loci_defence.tsv holds crude Fisher odds ratios.\n",
       "  Re-run 03_stats_loci.R (version 23/9) before drawing Figure 3.")

# =============================================================================
# ODDS RATIOS — each source parsed on its own terms, then combined.
# =============================================================================
MIN_CARRIERS <- 30

parse_or <- function(df, sec, klass) {
  d <- df %>% filter(section == sec)
  if (!nrow(d)) return(NULL)
  qq <- if ("q" %in% names(d)) d$q else
    suppressWarnings(as.numeric(sub(".*BH q = ", "", d$note)))
  tibble(Class = klass,
         category_final = sub(" ~ .*$", "", d$comparison),
         wing = case_when(grepl(" ~ Origin$", d$comparison) ~ "Origin",
                          grepl(" ~ Life$",   d$comparison) ~ "Life",
                          TRUE ~ NA_character_),
         OR = d$estimate, lo = d$ci_low, hi = d$ci_high, q = qq) %>%
    filter(!is.na(wing))
}
ORS <- bind_rows(parse_or(S_LOCI, "category enrichment (loci)", "Defence"),
                 parse_or(S_FIG3, "category enrichment (AVG)",  "AMG"))
stopifnot(nrow(ORS) > 0, any(ORS$Class == "AMG"), any(!is.na(ORS$q[ORS$Class == "AMG"])))
sig_cats <- ORS %>% filter(!is.na(q), q < 0.05) %>% pull(category_final) %>% unique()

# =============================================================================
# ROW ORDER — decided once, used by panels a and b
# =============================================================================
mb <- votu %>% group_by(Origin) %>%
  summarise(Mb = sum(contig_length) / 1e6, .groups = "drop")

keep <- long %>% distinct(Class, category_final, vOTU_ID) %>%
  count(Class, category_final, name = "carriers") %>%
  filter(carriers >= MIN_CARRIERS | category_final %in% sig_cats)

dens <- long %>% semi_join(keep, by = c("Class", "category_final")) %>%
  count(Class, category_final, Origin, Life, name = "n") %>%
  left_join(mb, by = "Origin") %>%
  mutate(per_Mb = n / Mb,
         Plot = ifelse(Origin == "Culture-independent", -per_Mb, per_Mb))

row_ord <- dens %>% group_by(Class, category_final, Origin) %>%
  summarise(pm = sum(per_Mb), .groups = "drop") %>%
  pivot_wider(names_from = Origin, values_from = pm, values_fill = 0) %>%
  mutate(ratio = `Culture-dependent` / pmax(`Culture-independent`, 1e-9)) %>%
  arrange(Class, ratio)
LEV <- row_ord$category_final
dens$category_final <- factor(dens$category_final, levels = LEV)

# =============================================================================
# PANEL A
# =============================================================================
xmax <- max(dens %>% group_by(Class, category_final, Origin) %>%
              summarise(t = sum(abs(Plot)), .groups = "drop") %>% pull(t))

p_a <- ggplot(dens, aes(category_final, Plot, fill = Life)) +
  geom_col(colour = "black", linewidth = 0.25, width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.7) +
  coord_flip(ylim = c(-xmax * 1.05, xmax * 1.05)) +
  facet_grid(Class ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = col_life, name = "Lifestyle") +
  scale_y_continuous(labels = abs, breaks = pretty_breaks(6),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = NULL,
       y = "per Mb of assembled viral sequence\n← culture-independent          culture-dependent →") +
  theme_pub +
  theme(strip.placement = "outside", legend.position = "none",
        axis.text.y = element_text(size = SZ_CATS))
life_legend <- grab_legend(p_a + theme(legend.position = "bottom"))

# =============================================================================
# PANEL B — Defence from the loci statistics, AMG from the gene statistics.
# Both are Mantel-Haenszel: each contrast with the other variable held
# constant, which is what the column headers now say.
# =============================================================================
fp <- ORS %>% filter(category_final %in% LEV) %>%
  mutate(q = ifelse(is.na(q), 1, q),
         category_final = factor(category_final, levels = LEV),
         Class = factor(Class, levels = c("Defence", "AMG")),
         Wing  = factor(ifelse(wing == "Origin",
                               "Culture-dep. vs -indep.\n(within lifestyles)",
                               "Temperate vs virulent\n(within fractions)"),
                        levels = c("Culture-dep. vs -indep.\n(within lifestyles)",
                                   "Temperate vs virulent\n(within fractions)")),
         lo = pmax(lo, 0.08), hi = pmin(hi, 30),
         sig = q < 0.05, lab = paste0(sprintf("%.1f", OR), star(q)))
stopifnot(nrow(fp) > 0)

p_b <- ggplot(fp, aes(OR, category_final)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.55) +
  geom_point(aes(fill = sig), shape = 21, size = 2.8, colour = "black") +
  geom_text(aes(x = 48, label = lab), hjust = 0, size = SZ_NUM) +
  facet_grid(Class ~ Wing, scales = "free_y", space = "free_y") +
  scale_x_log10(breaks = c(0.1, 1, 10), labels = c("0.1", "1", "10"),
                limits = c(0.07, 260)) +
  scale_fill_manual(values = c(`TRUE` = "#2C3E50", `FALSE` = "white"), guide = "none") +
  labs(x = "Mantel-Haenszel odds ratio (log scale)", y = NULL) +
  theme_pub +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        strip.text.y = element_blank(),
        panel.grid.major.y = element_line(colour = "grey94"))

# =============================================================================
# PANEL C — origin x lifestyle, in counts per Mb
# =============================================================================
inter <- votu %>% group_by(Origin, Life) %>%
  summarise(Mb = sum(contig_length) / 1e6,
            Defence = sum(n_Defence_loci), AMG = sum(n_AVG), .groups = "drop") %>%
  mutate(Defence = Defence / Mb, AMG = AMG / Mb) %>%
  pivot_longer(c(Defence, AMG), names_to = "Class", values_to = "value") %>%
  mutate(Class = factor(Class, levels = c("Defence", "AMG")))

lr <- bind_rows(
  S_LOCI %>% filter(section == "fraction vs lifestyle (loci)",
                    grepl("fraction beyond lifestyle", comparison)) %>%
    slice(1) %>% transmute(Class = "Defence", p),
  read.delim("stats_results.tsv", stringsAsFactors = FALSE) %>%
    filter(section == "fraction vs lifestyle", grepl("n_AVG", comparison)) %>%
    slice(1) %>% transmute(Class = "AMG", p)) %>%
  mutate(Class = factor(Class, levels = c("Defence", "AMG")),
         lab = paste0("origin, given lifestyle  ", ifelse(star(p) == "", "n.s.", star(p))))

p_c <- ggplot(inter, aes(Life, value, fill = Origin)) +
  geom_col(position = position_dodge(0.72), width = 0.62,
           colour = "black", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.1f", value)), position = position_dodge(0.72),
            vjust = -0.35, size = SZ_NUM) +
  geom_text(data = lr, aes(x = 1.5, y = Inf, label = lab), inherit.aes = FALSE,
            vjust = 1.5, size = SZ_SMALL, fontface = "italic", colour = "grey25") +
  facet_wrap(~Class, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = col_ori, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.30))) +
  labs(x = NULL, y = "per Mb") +
  theme_pub + theme(legend.position = "none")
origin_legend <- grab_legend(p_c + theme(legend.position = "bottom"))

# =============================================================================
# PANEL D — system families, faceted by axis with independent x scales
# =============================================================================
FAM_EN <- c("MTasa huerfana" = "Orphan methyltransferase",
            "Nucleasa putativa" = "Putative nuclease",
            "DNA-modification systems (other)" = "DNA modification (other)")
f3 <- loci %>%
  filter(completitud_locus == "completo", eje %in% c("1 Contra-defensa", "2 Defensa")) %>%
  mutate(Axis   = factor(ifelse(eje == "1 Contra-defensa", "Self-protection", "Exclusion"),
                         levels = c("Self-protection", "Exclusion")),
         Family = recode(familia, !!!FAM_EN))

TOP_PER_AXIS <- 6
fam <- f3 %>% count(Axis, Family, name = "n") %>%
  group_by(Axis) %>% arrange(desc(n), .by_group = TRUE) %>%
  mutate(Family = ifelse(row_number() <= TOP_PER_AXIS, Family,
                         paste0("Other (", n() - TOP_PER_AXIS, " families)"))) %>%
  count(Axis, Family, wt = n, name = "n") %>% ungroup() %>%
  mutate(Family = reorder(Family, n))

p_d <- ggplot(fam, aes(n, Family, fill = Axis)) +
  geom_col(colour = "black", linewidth = 0.3, width = 0.68) +
  geom_text(aes(label = n), hjust = -0.22, size = SZ_NUM, fontface = "bold") +
  facet_wrap(~Axis, ncol = 1, scales = "free") +
  scale_fill_manual(values = col_axis, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.26))) +
  labs(x = "complete loci", y = NULL) +
  theme_pub + theme(axis.text.y = element_text(size = SZ_CATS))

# =============================================================================
# PANEL E — axis composition across fractions
# =============================================================================
ax <- f3 %>%
  left_join(votu %>% select(vOTU_ID, Frac), by = c("contig_id" = "vOTU_ID")) %>%
  count(Frac, Axis) %>% group_by(Frac) %>%
  mutate(pct = 100 * n / sum(n), tot = sum(n)) %>% ungroup()

or_lab <- sapply(levels(votu$Frac), function(f) {
  if (f == "VIR") return("ref")
  s <- S_FIG3 %>% filter(section == "axis composition",
                         comparison == paste0(f, " vs VIR")) %>% slice(1)
  if (nrow(s) == 0) "" else sprintf("OR %.1f", s$estimate)
})
cl <- S_FIG3 %>% filter(section == "axis composition",
                        grepl("lifestyle, stratified by fraction", comparison)) %>% slice(1)

p_e <- ggplot(ax, aes(Frac, pct, fill = Axis)) +
  geom_col(colour = "black", linewidth = 0.3, width = 0.6) +
  geom_text(data = ax %>% filter(Axis == "Exclusion"),
            aes(y = 103, label = sprintf("%.1f%%", pct)), hjust = 0,
            fontface = "bold", size = SZ_NUM) +
  geom_text(data = ax %>% distinct(Frac, tot),
            aes(Frac, y = -2, label = paste0("n=", tot, "  ", or_lab[as.character(Frac)])),
            inherit.aes = FALSE, hjust = 1, size = SZ_SMALL, colour = "grey30") +
  annotate("text", x = 4.55, y = 50,
           label = paste0("lifestyle, given fraction  ",
                          ifelse(star(cl$p) == "", "n.s.", star(cl$p))),
           size = SZ_SMALL, fontface = "italic", colour = "grey25") +
  coord_flip(ylim = c(-40, 122), clip = "off") +
  scale_fill_manual(values = col_axis, name = NULL) +
  scale_y_continuous(breaks = seq(0, 100, 25), labels = function(x) paste0(x, "%")) +
  labs(x = "free particles  →  cultured cells", y = "% of complete loci") +
  theme_pub + theme(legend.position = "none")
axis_legend <- grab_legend(p_e + theme(legend.position = "bottom"))

# =============================================================================
# ASSEMBLY — letters drawn by cowplot so they sit at one height per row
# =============================================================================
LAB <- function(...) list(labels = c(...), label_size = BASE + 1,
                          label_fontface = "bold", label_x = 0.005, hjust = 0,
                          label_y = 0.995, vjust = 1)

row1 <- do.call(plot_grid, c(list(p_a, p_b, ncol = 2, rel_widths = c(1.45, 1),
                                  align = "h", axis = "tb"),
                             LAB("A | Density by category", "B | Enrichment")))
top  <- plot_grid(row1,
                  plot_grid(NULL, life_legend, NULL, ncol = 3,
                            rel_widths = c(0.35, 0.30, 0.35)),
                  ncol = 1, rel_heights = c(1, 0.07))

row2 <- do.call(plot_grid, c(list(p_c, p_d, p_e, ncol = 3,
                                  rel_widths = c(1.0, 1.0, 1.05),
                                  align = "h", axis = "tb"),
                             LAB("C | Origin and lifestyle combined",
                                 "D | System families of the complete loci",
                                 "E | Axis composition across fractions")))
bottom <- plot_grid(row2,
                    plot_grid(NULL, origin_legend, axis_legend, NULL, ncol = 4,
                              rel_widths = c(0.10, 0.40, 0.40, 0.10)),
                    ncol = 1, rel_heights = c(1, 0.09))

fig3 <- plot_grid(top, bottom, ncol = 1, rel_heights = c(1, 1))

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

export_isme(fig3, 3, 15, 11.5)

# =============================================================================
# REPORT
# =============================================================================
rule <- function(x) cat("\n", strrep("-", 74), "\n ", x, "\n", strrep("-", 74), "\n", sep = "")
rule("categories plotted");        print(as.data.frame(keep))
rule("a — density per Mb");        print(as.data.frame(dens %>% mutate(per_Mb = round(per_Mb, 2)) %>%
                                                         select(Class, category_final, Origin, Life, n, per_Mb) %>% arrange(Class, category_final)))
rule("b — odds ratios plotted");   print(as.data.frame(fp %>% select(Class, category_final, Wing, OR, lo, hi, q) %>%
                                                         mutate(across(c(OR, lo, hi), ~round(.x, 2)), q = signif(q, 3))))
rule("c — per Mb by origin and lifestyle")
print(as.data.frame(inter %>% mutate(value = round(value, 2))))
rule("d — ALL families, not only the plotted ones")
print(as.data.frame(f3 %>% count(Axis, Family, sort = TRUE)))
rule("e — axis by fraction");      print(as.data.frame(ax %>% mutate(pct = round(pct, 1))))
cat("\nFig3: the four files are in ", OUT, "\n", sep = "")
