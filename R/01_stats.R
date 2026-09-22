# =============================================================================
# 01_stats.R — EVERY statistic for the manuscript, in one place.
#
# Consolidates the analyses that used to live inside Fig2_V2.R and Fig3_V3.R.
# The figure scripts must not compute a single test: they read
# stats_results.tsv. In the old Fig3_V3.R the BH correction had to be
# duplicated inside the figure so the panel matched the text; if a test were
# added the vector order shifted and the figure silently showed the wrong
# p-value. That failure mode is now impossible.
#
# OUTPUTS
#   stats_report.txt    readable report; paste numbers from here
#   stats_results.tsv   one row per test -> read by 02_Fig_2.R and
#                       05_Fig_3.R (panel c). The tests reported in the paper
#                       are listed in Supplementary Tables S3 and S4.
#
# TWO CURRENCIES, both necessary
#   COUNTS per vOTU : whole catalogue, all four fractions, immune to MDA bias
#   CPM / TPM       : only CB22 and GLI were mapped (VIRFC and LIS are MDA-
#                     amplified), so abundance tests use 4 vs 4 samples and
#                     are UNDERPOWERED BY DESIGN. Report abundance as
#                     supporting evidence; the count-based tests on 11,100
#                     vOTUs carry the weight.
#
# CONFOUNDS, each with its own test
#   group size    -> everything per vOTU, never pooled gene counts
#   contig length -> negative binomial GLM with log(length); CMH by length bin
#   completeness  -> CMH by completeness quartile; high-quality-only Fisher
#   geography     -> Zone included in the count models and in the PERMANOVA
#   mitomycin     -> only LIS induced, so CB22 vs GLI is the clean contrast
#   MDA           -> balanced across groups by design; tested explicitly
#
# FIXED IN THIS VERSION (22/9)
#   * PERMANOVA p-values are exact (all 35 partitions) instead of 999 random
#     permutations, and the Zone + Origin model is reported for its R2 only.
#   * Relative paths, and the input tables are accepted with or without the
#     "_28_8" suffix: run from the folder that contains them.
#
# FIXED IN THE PREVIOUS VERSION
#   * adonis2(..., by = "terms"): without it, recent vegan tests the whole
#     model in one line, so "does origin add beyond zone" was never actually
#     answered and the values recorded were from the Residual row.
#   * Zone added to the main negative binomial models, because the LR test
#     shows it matters (p ~ 1e-5 for both gene classes).
#   * AVG/defence co-occurrence adjusted for contig length: both are likelier
#     on long contigs, so the raw odds ratio can be pure genome-size effect.
#   * per_cat() calls in section 13 (sensitivity) were missing the `klass`
#     argument entirely -- a leftover from before per_cat() took a class
#     parameter. That made R error immediately ("argument klass is missing,
#     with no default") instead of silently producing wrong output, which
#     is the one good thing about it, but it still needs fixing: both
#     sensitivity subsets now call per_cat() once for Defence and once for
#     AVG, exactly like section 10 does.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(MASS); library(vegan)
})
select <- dplyr::select

# Run from the folder that contains the input tables, e.g. setwd("path/to/R").

# PERMANOVA p-value by complete enumeration. With four libraries per group
# there are only choose(8, 4) / 2 = 35 distinct partitions, so the smallest
# attainable p-value is 1/35 = 0.029 and random permutations cannot resolve
# it (999 random permutations gave 0.024). Every partition is evaluated once
# with vegan's own pseudo-F, and p is the proportion with F >= observed.
exact_perm_p <- function(d, g) {
  g <- factor(g); n <- length(g); k <- sum(g == levels(g)[1])
  F_of <- function(gr) adonis2(d ~ gr, data = data.frame(gr = gr),
                               permutations = 0)$F[1]
  F_obs <- F_of(g)
  splits <- combn(n, k)
  splits <- splits[, splits[1, ] == 1, drop = FALSE]   # each partition once
  F_all <- apply(splits, 2, function(s) {
    gr <- rep(levels(g)[2], n); gr[s] <- levels(g)[1]
    F_of(factor(gr, levels = levels(g)))
  })
  c(p = mean(F_all >= F_obs - 1e-9), n_partitions = length(F_all))
}


# Input files. The scripts accept the names with or without the "_28_8" suffix,
# whichever exists in the working directory.
pick <- function(...) {
  f <- c(...); ok <- f[file.exists(f)]
  if (!length(ok)) stop("No file found: ", paste(f, collapse = " / "))
  cat("   input:", ok[1], "\n"); ok[1]
}
F_VOTU <- pick("07_votu_with_counts_28_8.tsv", "07_votu_with_counts.tsv")
F_LONG <- pick("07_long_table_FINAL_28_8.tsv", "07_long_table_FINAL.tsv")

votu_all <- read.delim(F_VOTU, stringsAsFactors = FALSE)
long_all <- read.delim(F_LONG, stringsAsFactors = FALSE)

# The new tables carry the WHOLE assembly plus a flag; the old ones came
# pre-filtered. Filtering happens here, once, and nowhere else.
votu <- votu_all %>% filter(Passes_catalog == "Yes")
long <- long_all %>% filter(Passes_catalog == "Yes",
                            gene_class %in% c("AVG", "Defence"))

# Only the 8 real sample columns. TPM_total and TPM_sum_* are derived and
# would otherwise be counted as extra samples.
SAMPLE_RE <- "^TPM_(DE|EL|MB|OR)_(CB22|GLI)$"
tpm_cols  <- grep(SAMPLE_RE, names(votu), value = TRUE)
stopifnot(length(tpm_cols) == 8)

votu <- votu %>%
  mutate(Life_Gen  = ifelse(grepl("^Temperate", Lifestyle),
                            "Temperate", "Virulent"),
         temperate = Life_Gen == "Temperate",
         lenbin    = cut(contig_length, c(0, 5e3, 1e4, 2e4, Inf),
                         labels = c("<5kb", "5-10kb", "10-20kb", ">=20kb")),
         Viral_Fam = str_extract(Taxonomy, "f_[^;]+"),
         Viral_Fam = ifelse(is.na(Viral_Fam), "Unclassified",
                            gsub("f_", "", Viral_Fam)),
         Host_Class = str_extract(Host_Prediction, "c__[^;]+"),
         Host_Class = ifelse(is.na(Host_Class) |
                               grepl("Unknown", Host_Prediction),
                             "Unknown", gsub("c__", "", Host_Class)))

RES <- list()
add <- function(section, test, comparison, n1 = NA, n2 = NA, estimate = NA,
                ci_low = NA, ci_high = NA, p = NA, q = NA, note = "") {
  # q is only supplied by tests that belong to a BH-corrected family
  # (per_cat()). Everything else stays NA and falls back to raw p when
  # the final `signif` column is built at the very end of this script --
  # see res_tbl below. Before this fix, q lived only as text inside
  # `note` ("AVG; q=1.4e-01"), so the figure's asterisks were computed
  # from raw p even for tests whose whole point was the BH correction:
  # e.g. AVG "Amino acid metabolism" has raw p = 0.044 (a lone "*") but
  # BH q = 0.144 (properly not significant once corrected for testing
  # 13 categories at once).
  RES[[length(RES) + 1]] <<- tibble(
    section = section, test = test, comparison = comparison, n1 = n1, n2 = n2,
    estimate = estimate, ci_low = ci_low, ci_high = ci_high, p = p, q = q,
    note = note)
}

sink("stats_report.txt", split = TRUE)
rule <- function(x) cat("\n", strrep("=", 78), "\n ", x, "\n",
                        strrep("=", 78), "\n", sep = "")

cat("Statistics — Antarctic virome manuscript\n")
cat("generated:", format(Sys.time()), "\n")
cat("catalogue:", nrow(votu), "vOTUs |", nrow(long), "AVG+defence genes\n")

fish <- function(d, var, lv_a, lv_b, outcome) {
  d <- d %>% filter(.data[[var]] %in% c(lv_a, lv_b))
  tb <- table(factor(d[[var]], levels = c(lv_b, lv_a)), d[[outcome]])
  # lv_a/lv_b travel with the table so show_fisher() can pull n1/n2 by NAME
  # instead of position -- see the fix note below.
  list(tb = tb, ft = fisher.test(tb), lv_a = lv_a, lv_b = lv_b)
}
show_fisher <- function(name, r, section, comparison, note = "") {
  cat("\n", name, "\n", sep = "")
  print(r$tb)
  cat(sprintf("   OR = %.3f   95%% CI %.3f-%.3f   p = %.3g\n",
              r$ft$estimate, r$ft$conf.int[1], r$ft$conf.int[2], r$ft$p.value))
  # FIXED: this used to read n1 = sum(r$tb[1, ]), n2 = sum(r$tb[2, ]) --
  # positional indexing into a table whose ROW ORDER is fixed by fish() as
  # c(lv_b, lv_a), i.e. row 1 is always the SECOND-named group. Every
  # comparison built with fish() therefore had n1/n2 silently swapped
  # relative to what the `comparison` label says (e.g. "Culture-dependent
  # vs Culture-independent" reported n1 = the Culture-INDEPENDENT count).
  # Caught by cross-checking Table S2's real Fraction/Origin counts against
  # this file's n1/n2 columns by hand -- every one of the 7 affected rows
  # (lifestyle x2, design effects x2, completeness-control HQ-only, gene
  # content x2) was reversed. The OR/CI/p themselves were never affected:
  # fisher.test() uses the whole 2x2 table regardless of how we later label
  # the margins. per_cat() below never had this bug because it already
  # indexes by row NAME (tb["Culture-dependent","TRUE"]) -- this makes
  # show_fisher() do the same, using the lv_a/lv_b names fish() now carries,
  # so it cannot silently flip again no matter how the table is ordered.
  add(section, "Fisher exact", comparison, n1 = sum(r$tb[r$lv_a, ]),
      n2 = sum(r$tb[r$lv_b, ]), estimate = unname(r$ft$estimate),
      ci_low = r$ft$conf.int[1], ci_high = r$ft$conf.int[2],
      p = r$ft$p.value, note = note)
}

# =============================================================================
rule("1. CATALOGUE DESCRIPTION")

by_frac <- votu %>%
  mutate(Total_TPM = rowSums(across(all_of(tpm_cols)))) %>%
  group_by(Fraction, Fraction_group, MDA, Mitomycin) %>%
  summarise(vOTUs = n(), TPM_sum = round(sum(Total_TPM), 1), .groups = "drop") %>%
  mutate(vOTUs_pct = round(100 * vOTUs / sum(vOTUs), 2),
         TPM_pct   = round(100 * TPM_sum / sum(TPM_sum), 2)) %>%
  arrange(desc(vOTUs))
print(as.data.frame(by_frac))

cat("\nBy zone:\n")
print(as.data.frame(votu %>% count(Zone, Fraction_group) %>%
                      pivot_wider(names_from = Fraction_group, values_from = n, values_fill = 0)))

cat("\nEvidence lines supporting each vOTU:\n")
print(as.data.frame(votu %>% count(N_evidence) %>%
                      mutate(pct = round(100 * n / sum(n), 1))))

cat("\nLength / completeness by fraction:\n")
print(as.data.frame(votu %>% group_by(Fraction) %>%
                      summarise(n = n(), median_len = median(contig_length),
                                median_completeness = round(median(completeness, na.rm = TRUE), 1),
                                .groups = "drop")))

# =============================================================================
rule("2. LIFESTYLE PARTITIONING")

r <- fish(votu, "Fraction_group", "Culture-dependent",
          "Culture-independent", "temperate")
show_fisher("Temperate ~ culture group (pooled; LIS is mitomycin-induced)",
            r, "lifestyle", "Culture-dependent vs Culture-independent")

r <- fish(votu, "Fraction", "GLI", "CB22", "temperate")
show_fisher("Temperate ~ GLI vs CB22 (both non-MDA, neither induced)",
            r, "lifestyle", "GLI vs CB22 (clean contrast)",
            note = "mitomycin-free and MDA-free")

cat("\nBy fraction, for the record:\n")
print(as.data.frame(votu %>% group_by(Fraction, MDA, Mitomycin) %>%
                      summarise(n = n(), temperate = sum(temperate),
                                pct = round(100 * mean(temperate), 1), .groups = "drop") %>%
                      arrange(desc(pct))))

cat("\nCochran-Mantel-Haenszel, stratified by contig length\n")
arr <- with(votu %>% filter(Fraction %in% c("GLI", "CB22")),
            table(factor(Fraction, levels = c("CB22", "GLI")), temperate, lenbin))
print(ftable(arr))
cmh <- mantelhaen.test(arr, exact = FALSE)
cat(sprintf("   common OR = %.3f   95%% CI %.3f-%.3f   p = %.3g\n",
            cmh$estimate, cmh$conf.int[1], cmh$conf.int[2], cmh$p.value))
add("lifestyle", "Cochran-Mantel-Haenszel", "GLI vs CB22, by length bin",
    estimate = unname(cmh$estimate), ci_low = cmh$conf.int[1],
    ci_high = cmh$conf.int[2], p = cmh$p.value, note = "length-adjusted")

# =============================================================================
rule("3. COMPLETENESS CONTROL — is 'temperate' just 'better assembled'?")
# Calling a phage temperate requires a detectable integrase, so a more
# complete genome is likelier to be called temperate. Length is handled
# above; CheckV completeness measures that worry more directly.

cmp <- votu %>% filter(!is.na(completeness))
print(as.data.frame(cmp %>% group_by(Life_Gen) %>%
                      summarise(n = n(), mean = round(mean(completeness), 1),
                                sd = round(sd(completeness), 1),
                                median = round(median(completeness), 1), .groups = "drop")))
w <- wilcox.test(completeness ~ Life_Gen, data = cmp)
cat(sprintf("   Wilcoxon p = %.3g\n", w$p.value))
add("completeness control", "Wilcoxon rank-sum",
    "completeness: temperate vs virulent", p = w$p.value)

sp <- suppressWarnings(cor.test(cmp$completeness,
                                as.integer(cmp$Life_Gen == "Temperate"),
                                method = "spearman"))
cat(sprintf("   Spearman rho = %.3f   p = %.3g\n", sp$estimate, sp$p.value))
add("completeness control", "Spearman correlation",
    "completeness vs temperate", estimate = unname(sp$estimate), p = sp$p.value)

cat("\nCMH stratified by completeness quartile (GLI vs CB22)\n")
cq <- cmp %>% filter(Fraction %in% c("GLI", "CB22")) %>%
  mutate(cbin = cut(completeness, quantile(completeness, 0:4 / 4),
                    include.lowest = TRUE))
arr2 <- with(cq, table(factor(Fraction, levels = c("CB22", "GLI")),
                       temperate, cbin))
print(ftable(arr2))
c2 <- mantelhaen.test(arr2, exact = FALSE)
cat(sprintf("   common OR = %.3f   95%% CI %.3f-%.3f   p = %.3g\n",
            c2$estimate, c2$conf.int[1], c2$conf.int[2], c2$p.value))
add("completeness control", "Cochran-Mantel-Haenszel",
    "GLI vs CB22, by completeness quartile", estimate = unname(c2$estimate),
    ci_low = c2$conf.int[1], ci_high = c2$conf.int[2], p = c2$p.value)

hq <- votu %>% filter(checkv_quality %in%
                        c("Medium-quality", "High-quality", "Complete"))
r <- fish(hq, "Fraction_group", "Culture-dependent",
          "Culture-independent", "temperate")
show_fisher("Temperate ~ culture group (medium-quality or better only)",
            r, "completeness control", "dependent vs independent, HQ only")
cat("\nCheckV quality x lifestyle:\n")
print(as.data.frame(votu %>% count(checkv_quality, Life_Gen) %>%
                      pivot_wider(names_from = Life_Gen, values_from = n, values_fill = 0)))

# =============================================================================
rule("4. WITHIN-GROUP CONTRASTS: does MDA or mitomycin drive anything?")
# The design is balanced for MDA (one MDA and one non-MDA fraction on each
# side) but not for mitomycin (LIS only). These contrasts isolate each.
# NOTE THE LOGIC: mitomycin INDUCES prophages, so LIS should look MORE
# temperate than GLI. If it looks LESS, whatever suppresses temperate
# detection there is strong enough to reverse the induction effect.

r <- fish(votu, "Fraction", "VIRFC", "CB22", "temperate")
show_fisher("MDA effect within culture-INDEPENDENT (VIRFC vs CB22)",
            r, "design effects", "VIRFC vs CB22 (MDA vs not)",
            note = "both culture-independent; isolates MDA")

r <- fish(votu, "Fraction", "LIS", "GLI", "temperate")
show_fisher("MDA + mitomycin within culture-DEPENDENT (LIS vs GLI)",
            r, "design effects", "LIS vs GLI (MDA+mitomycin vs neither)",
            note = "both culture-dependent; MDA and mitomycin confounded")

# =============================================================================
rule("5. ALPHA DIVERSITY (per sample; CB22 and GLI only, 4 vs 4)")

tpm_v <- votu %>%
  select(vOTU_ID, all_of(tpm_cols)) %>%
  pivot_longer(-vOTU_ID, names_to = "Sample", values_to = "TPM") %>%
  filter(TPM > 0) %>%
  mutate(Zone   = sapply(strsplit(Sample, "_"), `[`, 2),
         FracS  = sapply(strsplit(Sample, "_"), `[`, 3),
         Origin = ifelse(FracS == "CB22",
                         "Culture-independent", "Culture-dependent"))

alpha <- tpm_v %>% group_by(Sample, Origin) %>%
  summarise(Richness = n(), Shannon = diversity(TPM, index = "shannon"),
            .groups = "drop")
print(as.data.frame(alpha %>% group_by(Origin) %>%
                      summarise(n = n(), Richness_mean = round(mean(Richness), 1),
                                Richness_sd = round(sd(Richness), 1),
                                Shannon_mean = round(mean(Shannon), 2),
                                Shannon_sd = round(sd(Shannon), 2), .groups = "drop")))
for (v in c("Richness", "Shannon")) {
  w <- wilcox.test(alpha[[v]] ~ alpha$Origin)
  cat(sprintf("   %-9s Wilcoxon p = %.4g\n", v, w$p.value))
  add("alpha diversity", "Wilcoxon rank-sum",
      paste0(v, ": dependent vs independent"), p = w$p.value,
      note = "4 vs 4 samples; minimum attainable p ~ 0.029")
}

# =============================================================================
rule("6. BETA DIVERSITY (community structure)")
# p-values are exact (exact_perm_p, 35 partitions). The design is balanced
# (4 stations x 2 origins), so the R2 of origin is identical with or without
# station in the model. The station model is printed for its R2 only and is
# not tested: permutations restricted within stations have a floor of 0.125.

otu <- tpm_v %>% select(Sample, vOTU_ID, TPM) %>%
  pivot_wider(names_from = vOTU_ID, values_from = TPM, values_fill = 0) %>%
  column_to_rownames("Sample")
meta <- tibble(Sample = rownames(otu)) %>%
  mutate(Origin = ifelse(grepl("_CB22$", Sample),
                         "Culture-independent", "Culture-dependent"),
         Zone = sapply(strsplit(Sample, "_"), `[`, 2)) %>% as.data.frame()
rownames(meta) <- meta$Sample

for (idx in c("bray", "jaccard")) {
  d <- if (idx == "jaccard") vegdist(ifelse(otu > 0, 1, 0), "jaccard")
  else vegdist(otu, "bray")
  cat("\n---", toupper(idx), "---\n")
  
  a1 <- adonis2(d ~ Origin, data = meta, permutations = 0, by = "terms")
  ex <- exact_perm_p(d, meta$Origin)
  print(a1)
  cat(sprintf("   Origin: R2 = %.3f  F = %.3f  exact p = %.4f (%d partitions)\n",
              a1$R2[1], a1$F[1], ex["p"], ex["n_partitions"]))
  add("beta diversity", paste0("PERMANOVA (", idx, ")"),
      "community ~ origin", estimate = a1$R2[1], p = ex["p"],
      note = "exact p over all 35 partitions, 8 samples")

  a2 <- adonis2(d ~ Zone + Origin, data = meta, permutations = 0, by = "terms")
  cat(sprintf("   R2 of origin with station in the model: %.3f (not tested)\n",
              a2$R2[rownames(a2) == "Origin"]))
  
  bd <- betadisper(d, meta$Origin)
  set.seed(1); pt <- permutest(bd, permutations = 999)
  cat("   betadisper p =", signif(pt$tab$`Pr(>F)`[1], 4),
      "  (NON-significant is what you want: PERMANOVA is then detecting a\n",
      "   location difference, not unequal dispersion)\n")
  add("beta diversity", paste0("betadisper (", idx, ")"),
      "homogeneity of dispersion", p = pt$tab$`Pr(>F)`[1])
}

# =============================================================================
rule("7. DARK MATTER AND BIOGEOGRAPHY")

dark <- tpm_v %>% left_join(votu %>% select(vOTU_ID, Viral_Fam, Host_Class),
                            by = "vOTU_ID") %>%
  group_by(Origin) %>%
  summarise(TPM = round(sum(TPM), 1),
            Viral_unclassified_pct =
              round(100 * sum(TPM[Viral_Fam == "Unclassified"]) / sum(TPM), 2),
            Host_unknown_pct =
              round(100 * sum(TPM[Host_Class == "Unknown"]) / sum(TPM), 2),
            .groups = "drop")
print(as.data.frame(dark))
cat("\nCatalogue level (not abundance-weighted):\n")
cat(sprintf("   unclassified family: %d (%.1f%%)\n",
            sum(votu$Viral_Fam == "Unclassified"),
            100 * mean(votu$Viral_Fam == "Unclassified")))
cat(sprintf("   unknown host       : %d (%.1f%%)\n",
            sum(votu$Host_Class == "Unknown"),
            100 * mean(votu$Host_Class == "Unknown")))

nz <- votu %>%
  mutate(n_zones = rowSums(sapply(c("DE", "EL", "MB", "OR"), function(z)
    rowSums(across(all_of(grep(paste0("^TPM_", z, "_"), tpm_cols,
                               value = TRUE)))) > 0))) %>%
  filter(n_zones > 0)
cat("\nZones per detected vOTU:\n")
print(as.data.frame(nz %>% count(n_zones) %>%
                      mutate(pct = round(100 * n / sum(n), 2))))

# =============================================================================
rule("8. GENE CONTENT — COUNTS (all four fractions, MDA-independent)")

print(as.data.frame(votu %>% group_by(Fraction_group) %>%
                      summarise(vOTUs = n(), defence_per_votu = round(mean(n_Defence), 3),
                                avg_per_votu = round(mean(n_AVG), 3),
                                pct_with_defence = round(100 * mean(n_Defence > 0), 1),
                                pct_with_avg = round(100 * mean(n_AVG > 0), 1), .groups = "drop")))

cat("\nNormalised per 100 vOTUs, by fraction:\n")
print(as.data.frame(votu %>% group_by(Fraction) %>%
                      summarise(vOTUs = n(),
                                defence_per_100 = round(100 * sum(n_Defence) / n(), 1),
                                avg_per_100 = round(100 * sum(n_AVG) / n(), 1), .groups = "drop")))

for (v in c("n_Defence", "n_AVG")) {
  d <- votu %>% mutate(has = .data[[v]] > 0)
  r <- fish(d, "Fraction_group", "Culture-dependent",
            "Culture-independent", "has")
  show_fisher(paste0(v, " > 0 ~ culture group"), r, "gene content",
              paste0(v, " presence: dependent vs independent"))
}

cat("\n--- negative binomial GLM, adjusted for length AND zone ---\n")
cat("    Zone is included because the LR test below shows it matters.\n")
for (v in c("n_Defence", "n_AVG")) {
  m <- try(glm.nb(as.formula(paste(v,
                                   "~ Fraction_group + Zone + log(contig_length)")), data = votu),
           silent = TRUE)
  if (inherits(m, "try-error")) next
  s <- summary(m)$coefficients
  ro <- grep("Fraction_group", rownames(s), value = TRUE)[1]
  cat("\n  ", v, " (reference = Culture-dependent, zone DE)\n", sep = "")
  print(round(s, 4))
  cat(sprintf("   IRR for culture-independent = %.3f   p = %.3g\n",
              exp(s[ro, 1]), s[ro, 4]))
  add("gene content", "Negative binomial GLM",
      paste0(v, " ~ Fraction_group + Zone + log(length)"),
      estimate = exp(s[ro, 1]), p = s[ro, 4],
      note = "IRR, adjusted for length and zone")
  
  m0 <- glm.nb(as.formula(paste(v, "~ Fraction_group + log(contig_length)")),
               data = votu)
  an <- anova(m0, m)
  cat("   adding Zone: LR p =", signif(an$`Pr(Chi)`[2], 4), "\n")
  add("gene content", "Likelihood-ratio test",
      paste0(v, ": Zone beyond fraction and length"), p = an$`Pr(Chi)`[2])
}

# =============================================================================
rule("9. IS IT THE FRACTION, OR JUST BEING TEMPERATE?")

for (v in c("n_Defence", "n_AVG")) {
  m <- try(glm.nb(as.formula(paste(v,
                                   "~ Fraction_group * Life_Gen + Zone + log(contig_length)")),
                  data = votu), silent = TRUE)
  if (inherits(m, "try-error")) next
  cat("\n  ", v, "\n", sep = ""); print(round(summary(m)$coefficients, 4))
  m0 <- glm.nb(as.formula(paste(v, "~ Life_Gen + Zone + log(contig_length)")),
               data = votu)
  an <- anova(m0, m)
  cat("   adding Fraction_group to a model that already has Life_Gen:",
      "LR p =", signif(an$`Pr(Chi)`[2], 4), "\n")
  add("fraction vs lifestyle", "Likelihood-ratio test",
      paste0(v, ": fraction beyond lifestyle"), p = an$`Pr(Chi)`[2],
      note = "NB GLM, adjusted for length and zone")
}
cat("\nCross-tabulation:\n")
print(as.data.frame(votu %>% group_by(Fraction_group, Life_Gen) %>%
                      summarise(vOTUs = n(), defence_per_votu = round(mean(n_Defence), 3),
                                avg_per_votu = round(mean(n_AVG), 3), .groups = "drop")))

# =============================================================================
rule("10. PER-CATEGORY ENRICHMENT (each its own 2x2 + BH)")

per_cat <- function(long_df, votu_df, label, section, klass) {
  carriers <- long_df %>% filter(gene_class == klass) %>%
    distinct(votu_id, category_final)
  base <- votu_df %>% select(vOTU_ID, Fraction_group)
  res <- map_dfr(sort(unique(carriers$category_final)), function(cc) {
    has <- base %>% mutate(has = vOTU_ID %in%
                             carriers$votu_id[carriers$category_final == cc])
    tb <- table(factor(has$Fraction_group,
                       levels = c("Culture-independent", "Culture-dependent")),
                has$has)
    if (any(dim(tb) < 2)) return(NULL)
    ft <- fisher.test(tb)
    tibble(gene_class = klass, category = cc,
           n_dep = tb["Culture-dependent", "TRUE"],
           n_indep = tb["Culture-independent", "TRUE"],
           pct_dep = round(100 * tb["Culture-dependent", "TRUE"] /
                             sum(tb["Culture-dependent", ]), 2),
           pct_indep = round(100 * tb["Culture-independent", "TRUE"] /
                               sum(tb["Culture-independent", ]), 2),
           OR = unname(ft$estimate),
           lo = ft$conf.int[1], hi = ft$conf.int[2],   # <- FALTABAN
           p = ft$p.value)
  }) %>% mutate(q = p.adjust(p, method = "BH")) %>% arrange(desc(OR))
  cat("\n", label, " [", klass, "]\n", sep = "")
  print(as.data.frame(res %>% mutate(across(where(is.numeric), ~signif(.x, 4)))))
  for (i in seq_len(nrow(res)))
    add(section, "Fisher exact + BH", res$category[i], n1 = res$n_dep[i],
        n2 = res$n_indep[i], estimate = res$OR[i],
        ci_low = res$lo[i], ci_high = res$hi[i], p = res$p[i], q = res$q[i],
        note = paste0(klass, "; q=", signif(res$q[i], 4)))
  res
}

# BH is applied within each gene class: the two families of tests answer
# different questions and should not dilute each other.
invisible(per_cat(long, votu, "ALL GENES", "category enrichment", "Defence"))
invisible(per_cat(long, votu, "ALL GENES", "category enrichment", "AVG"))

cat("\nAVG categories by culture group:\n")
print(as.data.frame(long %>% filter(gene_class == "AVG") %>%
                      count(category_final, Fraction_group) %>%
                      pivot_wider(names_from = Fraction_group, values_from = n, values_fill = 0)))

# =============================================================================
rule("11. DO AVGs AND DEFENCE GENES CO-OCCUR?")
# Both are likelier on long contigs, so a raw odds ratio can be pure
# genome-size effect. The logistic model below is the one to report.

co <- votu %>% mutate(hasA = n_AVG > 0, hasD = n_Defence > 0)
tb <- table(co$hasA, co$hasD, dnn = c("has AVG", "has defence"))
print(tb)
ft <- fisher.test(tb)
cat(sprintf("   crude OR = %.3f   95%% CI %.3f-%.3f   p = %.3g\n",
            ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))
add("co-occurrence", "Fisher exact (crude)",
    "AVG and defence in the same vOTU", estimate = unname(ft$estimate),
    ci_low = ft$conf.int[1], ci_high = ft$conf.int[2], p = ft$p.value,
    note = "NOT length-adjusted; see the model below")

mco <- glm(hasD ~ hasA + log(contig_length) + Fraction_group,
           data = co, family = binomial)
s <- summary(mco)$coefficients
print(round(s, 4))
ci <- suppressMessages(confint(mco))
cat(sprintf("   ADJUSTED OR for hasA = %.3f   95%% CI %.3f-%.3f   p = %.3g\n",
            exp(s["hasATRUE", 1]), exp(ci["hasATRUE", 1]),
            exp(ci["hasATRUE", 2]), s["hasATRUE", 4]))
cat("   If the adjusted OR falls towards 1, the co-occurrence was a\n")
cat("   genome-size artefact and the two strategies are independent.\n")
cat("   If it stays high, some phages genuinely carry both.\n")
add("co-occurrence", "Logistic regression (adjusted)",
    "defence ~ AVG + log(length) + fraction",
    estimate = exp(s["hasATRUE", 1]), ci_low = exp(ci["hasATRUE", 1]),
    ci_high = exp(ci["hasATRUE", 2]), p = s["hasATRUE", 4],
    note = "adjusted for contig length and fraction")

# =============================================================================
rule("12. ABUNDANCE (CPM) — CB22 and GLI only, 4 vs 4, UNDERPOWERED")

tpm_g <- long %>%
  pivot_longer(all_of(intersect(tpm_cols, names(long))),
               names_to = "Sample", values_to = "TPM") %>%
  filter(TPM > 0) %>%
  mutate(FracS  = sapply(strsplit(Sample, "_"), `[`, 3),
         Origin = ifelse(FracS == "CB22",
                         "Culture-independent", "Culture-dependent"))

by_sample <- tpm_g %>% group_by(Sample, Origin, gene_class) %>%
  summarise(TPM_sum = sum(TPM), .groups = "drop")
print(as.data.frame(by_sample %>% group_by(Origin, gene_class) %>%
                      summarise(n_samples = n(), mean_TPM = round(mean(TPM_sum), 1),
                                median_TPM = round(median(TPM_sum), 1), .groups = "drop")))

for (gt in c("AVG", "Defence")) {
  s <- by_sample %>% filter(gene_class == gt)
  if (n_distinct(s$Origin) < 2) next
  w <- wilcox.test(TPM_sum ~ Origin, data = s)
  cat(sprintf("   %-8s Wilcoxon p = %.4g\n", gt, w$p.value))
  add("abundance", "Wilcoxon rank-sum",
      paste0("CPM of ", gt, ": dependent vs independent"), p = w$p.value,
      note = "4 vs 4 samples; minimum attainable p ~ 0.029")
}

cat("\nCPM by category:\n")
print(as.data.frame(tpm_g %>% group_by(gene_class, category_final) %>%
                      summarise(TPM_sum = round(sum(TPM), 1), n = n(), .groups = "drop") %>%
                      arrange(gene_class, desc(TPM_sum))))

cat("\n--- PERMANOVA on functional composition ---\n")
mat <- tpm_g %>% group_by(Sample, category_final) %>%
  summarise(TPM = sum(TPM), .groups = "drop") %>%
  pivot_wider(names_from = category_final, values_from = TPM, values_fill = 0) %>%
  column_to_rownames("Sample")
mf <- tibble(Sample = rownames(mat)) %>%
  mutate(Origin = ifelse(grepl("_CB22$", Sample),
                         "Culture-independent", "Culture-dependent")) %>%
  as.data.frame()
d_fun <- vegdist(mat, "bray")
adf <- adonis2(d_fun ~ Origin, data = mf, permutations = 0, by = "terms")
ex  <- exact_perm_p(d_fun, mf$Origin)
print(adf)
cat(sprintf("   Origin: R2 = %.3f  exact p = %.4f\n", adf$R2[1], ex["p"]))
add("abundance", "PERMANOVA (Bray-Curtis)", "functional composition ~ origin",
    estimate = adf$R2[1], p = ex["p"], note = "exact p over all 35 partitions, 8 samples")

# =============================================================================
rule("13. SENSITIVITY — the conclusion must survive both restrictions")

print(as.data.frame(bind_rows(
  long %>% count(gene_class) %>% mutate(subset = "all genes"),
  long %>% filter(passes_context == 1) %>% count(gene_class) %>%
    mutate(subset = "aux_score <= 3"),
  long %>% filter(Conservative_10kb == "Yes") %>% count(gene_class) %>%
    mutate(subset = ">= 10 kb or circular")) %>%
    pivot_wider(names_from = gene_class, values_from = n)))

for (klass in c("Defence", "AVG")) {
  invisible(per_cat(long %>% filter(passes_context == 1), votu,
                    "RESTRICTED TO aux_score <= 3", "sensitivity aux<=3", klass))
  invisible(per_cat(long %>% filter(Conservative_10kb == "Yes"),
                    votu %>% filter(Conservative_10kb == "Yes"),
                    "RESTRICTED TO >= 10 kb OR CIRCULAR", "sensitivity 10kb", klass))
}

# =============================================================================
rule("14. NORMALISATION — per 100 vOTUs vs per Mb of assembled sequence")
# Genes scale almost linearly with contig length (log(length) coefficient
# ~0.95), and culture-dependent vOTUs are longer. So "genes per 100 vOTUs"
# still favours the group with bigger genomes. Genes per Mb removes both
# effects and is the descriptive equivalent of the length-adjusted GLM.

norm <- votu %>% group_by(Fraction_group) %>%
  summarise(vOTUs = n(),
            Mb = sum(contig_length) / 1e6,
            median_len = median(contig_length),
            defence = sum(n_Defence), avg = sum(n_AVG), .groups = "drop") %>%
  mutate(defence_per_100_votu = round(100 * defence / vOTUs, 1),
         avg_per_100_votu     = round(100 * avg / vOTUs, 1),
         defence_per_Mb       = round(defence / Mb, 2),
         avg_per_Mb           = round(avg / Mb, 2))
print(as.data.frame(norm))

cat("\nratio dependent / independent under each denominator:\n")
r100d <- norm$defence_per_100_votu[norm$Fraction_group == "Culture-dependent"] /
  norm$defence_per_100_votu[norm$Fraction_group == "Culture-independent"]
rMbd  <- norm$defence_per_Mb[norm$Fraction_group == "Culture-dependent"] /
  norm$defence_per_Mb[norm$Fraction_group == "Culture-independent"]
r100a <- norm$avg_per_100_votu[norm$Fraction_group == "Culture-dependent"] /
  norm$avg_per_100_votu[norm$Fraction_group == "Culture-independent"]
rMba  <- norm$avg_per_Mb[norm$Fraction_group == "Culture-dependent"] /
  norm$avg_per_Mb[norm$Fraction_group == "Culture-independent"]
cat(sprintf("   defence  per 100 vOTUs %.2fx   per Mb %.2fx\n", r100d, rMbd))
cat(sprintf("   AVG      per 100 vOTUs %.2fx   per Mb %.2fx\n", r100a, rMba))
cat("   The per-Mb ratio is the honest one. If it is much smaller than the\n")
cat("   per-vOTU ratio, part of the enrichment was genome size.\n")

cat("\nBy category, both denominators:\n")
mb <- votu %>% group_by(Fraction_group) %>%
  summarise(Mb = sum(contig_length) / 1e6, vOTUs = n(), .groups = "drop")
print(as.data.frame(
  long %>% filter(gene_class == "Defence") %>%
    count(category_final, Fraction_group) %>%
    left_join(mb, by = "Fraction_group") %>%
    mutate(per_100_votu = round(100 * n / vOTUs, 2),
           per_Mb = round(n / Mb, 2)) %>%
    select(category_final, Fraction_group, n, per_100_votu, per_Mb) %>%
    arrange(category_final, Fraction_group)))


# =============================================================================
rule("DONE")
res_tbl <- bind_rows(RES) %>%
  mutate(sig_basis = coalesce(q, p),   # BH-corrected q where it exists
         # (category enrichment, both
         # sensitivity subsets); raw p for
         # every other test, which was
         # never part of a multiple-testing
         # family in the first place
         across(where(is.numeric), ~signif(.x, 5)),
         signif = case_when(is.na(sig_basis) ~ "", sig_basis < 0.001 ~ "***",
                            sig_basis < 0.01 ~ "**", sig_basis < 0.05 ~ "*",
                            TRUE ~ "ns")) %>%
  select(-sig_basis)
write_tsv(res_tbl, "stats_results.tsv")
cat("\n", nrow(res_tbl), " tests written to stats_results.tsv\n", sep = "")
cat("Figure scripts read that file and never recompute a test.\n")
sink()
cat("\nWritten to", getwd(), "\n")
