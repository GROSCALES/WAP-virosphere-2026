# =============================================================================
# 04_stats_fig3.R — every statistic behind Figure 3 that is not the defence
# loci enrichment (that one is in 03_stats_loci.R).
#
# VERSION 23/9. THREE CHANGES over the 28/8 script:
#
#   1. THE PER-CATEGORY ENRICHMENT IS MANTEL-HAENSZEL. Same reasoning as in
#      03_stats_loci.R: origin and lifestyle are confounded, so each contrast
#      is estimated with the other held constant. AMG categories are counted
#      on genes and defence categories on genes here (the figure uses the
#      loci version from 03_stats_loci.R for defence, and this file only for
#      the AMG rows). The crude Fisher values are still written, under
#      "category enrichment (AVG, unadjusted)" / "(Defence, unadjusted)".
#
#   2. Integer-overflow guard: mantelhaen.test does its arithmetic in
#      integers and returned NA confidence intervals on tables of 11,100
#      vOTUs. Tables are cast to double before the test, which is why the
#      carrier-prevalence rows now carry the intervals printed in Table S4c.
#
#   3. No absolute path; input names tried with and without "_28_8".
#
# The three questions, which have three different answers:
#   Q1  how much defence does a vOTU carry?   -> lifestyle YES, origin YES
#   Q2  how many AMGs does it carry?          -> lifestyle NO,  origin YES
#   Q3  what KIND of defence is it (axis)?    -> lifestyle NO,  fraction YES
#
# OUTPUTS  stats_results_fig3.tsv, fig3_panel_data.tsv, fig3_butterfly.tsv,
#          fig3_axis_ladder.tsv, fig3_families.tsv, fig3_demography.tsv,
#          stats_report_fig3.txt
# =============================================================================

suppressPackageStartupMessages({library(tidyverse)})
select <- dplyr::select

pick <- function(...) {
  f <- c(...); ok <- f[file.exists(f)]
  if (!length(ok)) stop("No file found: ", paste(f, collapse = " / "))
  cat("   input:", ok[1], "\n"); ok[1]
}
votu <- read.delim(pick("07_votu_with_counts_28_8.tsv", "07_votu_with_counts.tsv"),
                   stringsAsFactors = FALSE) %>% filter(Passes_catalog == "Yes")
long <- read.delim(pick("07_long_table_FINAL_28_8.tsv", "07_long_table_FINAL.tsv"),
                   stringsAsFactors = FALSE) %>% filter(Passes_catalog == "Yes")
loci <- read.delim(pick("10_loci_defensa_28_8.tsv", "10_loci_defensa.tsv"),
                   stringsAsFactors = FALSE)

tpm_cols <- grep("^TPM_(DE|EL|MB|OR)_(CB22|GLI)$", names(votu), value = TRUE)
stopifnot(length(tpm_cols) == 8)
FR <- c(CB22 = "BAC", GLI = "CUL", VIRFC = "VIR", LIS = "LIS")
FR_LEV <- c("VIR", "LIS", "BAC", "CUL")

votu <- votu %>%
  mutate(Life   = factor(ifelse(grepl("^Temperate", Lifestyle), "Temperate", "Virulent"),
                         levels = c("Virulent", "Temperate")),
         Frac   = factor(FR[Fraction], levels = FR_LEV),
         Origin = factor(ifelse(Fraction_group == "Culture-dependent", "CD", "CI"),
                         levels = c("CI", "CD")),
         loglen = log10(contig_length),
         Total_TPM = rowSums(across(all_of(tpm_cols))))
long <- long %>%
  left_join(votu %>% select(vOTU_ID, Life, Frac, Origin), by = "vOTU_ID",
            suffix = c("", ".v"))

sink("stats_report_fig3.txt", split = TRUE)
rule <- function(x) cat("\n", strrep("=", 78), "\n ", x, "\n", strrep("=", 78), "\n", sep = "")
cat("Figure 3 statistics —", format(Sys.time()), "\n")
cat(nrow(votu), "vOTUs |", nrow(long), "classified genes |", nrow(loci), "loci\n")

RES <- list()
add <- function(section, test, comparison, n1 = NA, n2 = NA, estimate = NA,
                ci_low = NA, ci_high = NA, p = NA, note = "") {
  RES[[length(RES) + 1]] <<- tibble(section, test, comparison, n1, n2,
                                    estimate, ci_low, ci_high, p, note)
}
# as.double: see header, change 2
as_num_table <- function(tb) array(as.double(tb), dim(tb), dimnames(tb))
cmh_report <- function(arr, label, section) {
  r <- mantelhaen.test(as_num_table(arr), exact = FALSE)
  cat(sprintf("   %-46s OR = %6.3f  95%% CI %.2f-%.2f  p = %.4g\n",
              label, r$estimate, r$conf.int[1], r$conf.int[2], r$p.value))
  add(section, "Cochran-Mantel-Haenszel", label, estimate = unname(r$estimate),
      ci_low = r$conf.int[1], ci_high = r$conf.int[2], p = r$p.value)
  invisible(r)
}

# =============================================================================
rule("1. CARRIER PREVALENCE — does a vOTU carry any defence gene / any AMG?")
prev <- votu %>% group_by(Frac, Life) %>%
  summarise(n_votus = n(), Mb = sum(contig_length) / 1e6,
            defence_n = sum(n_Defence), amg_n = sum(n_AVG),
            pct_votus_def = 100 * mean(n_Defence > 0),
            pct_votus_amg = 100 * mean(n_AVG > 0), .groups = "drop") %>%
  mutate(Defence_per_Mb = defence_n / Mb, AMG_per_Mb = amg_n / Mb)
print(as.data.frame(prev %>% mutate(across(where(is.numeric), ~round(.x, 2)))))

cat("\nIs the lifestyle effect real, or the fraction in disguise?\n")
for (v in c("n_Defence", "n_AVG")) {
  cat("\n  ", v, "\n", sep = "")
  d <- votu %>% mutate(has = .data[[v]] > 0)
  klass <- ifelse(v == "n_Defence", "defence", "AMG")
  cmh_report(with(d, table(Life, has, Frac)),
             paste0("carries ", klass, " ~ lifestyle, stratified by fraction"),
             "carrier prevalence")
  cmh_report(with(d, table(Origin, has, Life)),
             paste0("carries ", klass, " ~ culture-dep vs -indep, stratified by lifestyle"),
             "carrier prevalence")
  cat("     stratum ORs for lifestyle:\n")
  for (f in FR_LEV) {
    s <- d %>% filter(Frac == f); ft <- fisher.test(table(s$Life, s$has))
    cat(sprintf("       %-4s T %5.1f%%  V %5.1f%%   OR = %5.2f  p = %.3g\n", f,
                100 * mean(s$has[s$Life == "Temperate"]),
                100 * mean(s$has[s$Life == "Virulent"]), ft$estimate, ft$p.value))
  }
}

# =============================================================================
rule("2. PER-CATEGORY CARRIER ENRICHMENT (Fig. 3b) — Mantel-Haenszel")
carrier_or <- function(klass, var, strat, lv_lo, lv_hi) {
  carr <- long %>% filter(gene_class == klass) %>% distinct(vOTU_ID, category_final)
  map_dfr(sort(unique(carr$category_final)), function(cc) {
    d <- votu %>% mutate(has = vOTU_ID %in% carr$vOTU_ID[carr$category_final == cc])
    tb <- as_num_table(table(factor(d[[var]], levels = c(lv_lo, lv_hi)), d$has, d[[strat]]))
    r  <- try(mantelhaen.test(tb, exact = FALSE), silent = TRUE)
    ft <- fisher.test(table(factor(d[[var]], levels = c(lv_lo, lv_hi)), d$has))
    ok <- !inherits(r, "try-error") && is.finite(r$estimate)
    tibble(gene_class = klass, split = var, category = cc,
           n_a = sum(d$has & d[[var]] == lv_hi), n_b = sum(d$has & d[[var]] == lv_lo),
           OR = if (ok) unname(r$estimate) else unname(ft$estimate),
           lo = if (ok) r$conf.int[1] else ft$conf.int[1],
           hi = if (ok) r$conf.int[2] else ft$conf.int[2],
           p  = if (ok) r$p.value else ft$p.value,
           OR_crude = unname(ft$estimate), p_crude = ft$p.value,
           test = if (ok) "CMH" else "Fisher (CMH undefined: empty stratum)")
  }) %>% mutate(q = p.adjust(p, "BH"), q_crude = p.adjust(p_crude, "BH")) %>%
    arrange(desc(n_a + n_b))
}
BUT <- bind_rows(
  carrier_or("Defence", "Origin", "Life", "CI", "CD"),
  carrier_or("Defence", "Life",   "Frac", "Virulent", "Temperate"),
  carrier_or("AVG",     "Origin", "Life", "CI", "CD"),
  carrier_or("AVG",     "Life",   "Frac", "Virulent", "Temperate"))
for (k in c("Defence", "AVG")) for (s in c("Origin", "Life")) {
  cat("\n---", k, "by", s, "(Mantel-Haenszel; crude Fisher alongside) ---\n")
  print(as.data.frame(BUT %>% filter(gene_class == k, split == s) %>%
                        mutate(across(c(OR, lo, hi, OR_crude), ~round(.x, 3)),
                               across(c(p, q, p_crude), ~signif(.x, 3)))))
}
for (i in seq_len(nrow(BUT))) with(BUT[i, ], {
  add(paste0("category enrichment (", gene_class, ")"), test,
      paste0(category, " ~ ", split), n1 = n_a, n2 = n_b, estimate = OR,
      ci_low = lo, ci_high = hi, p = p, note = paste0("BH q = ", signif(q, 3)))
  add(paste0("category enrichment (", gene_class, ", unadjusted)"), "Fisher exact",
      paste0(category, " ~ ", split), n1 = n_a, n2 = n_b, estimate = OR_crude,
      p = p_crude, note = paste0("crude; BH q = ", signif(q_crude, 3)))
})

# =============================================================================
rule("3. ABUNDANCE CARRIED — descriptive, 8 libraries only")
cpm <- votu %>% group_by(Frac, Life) %>%
  summarise(CPM_total = sum(Total_TPM),
            CPM_defence = sum(Total_TPM[n_Defence > 0]),
            CPM_amg     = sum(Total_TPM[n_AVG > 0]), .groups = "drop") %>%
  mutate(pct_defence = 100 * CPM_defence / CPM_total,
         pct_amg     = 100 * CPM_amg / CPM_total)
print(as.data.frame(cpm %>% mutate(across(where(is.numeric), ~round(.x, 1)))))

# =============================================================================
rule("4. WHAT KIND OF DEFENCE — axis composition (Fig. 3d, e)")
f3 <- loci %>%
  filter(completitud_locus == "completo", eje %in% c("1 Contra-defensa", "2 Defensa")) %>%
  left_join(votu %>% select(vOTU_ID, Frac, Life, Origin, loglen),
            by = c("contig_id" = "vOTU_ID")) %>%
  mutate(ax2 = as.integer(eje == "2 Defensa"))
stopifnot(!any(is.na(f3$Frac)))
cat("complete loci on axes 1-2:", nrow(f3), "\n")
print(as.data.frame(loci %>% count(completitud_locus)))
print(as.data.frame(f3 %>% count(eje, metrica)))
cat("\nfamilies per axis:", n_distinct(f3$familia[f3$ax2 == 0]), "self-protection,",
    n_distinct(f3$familia[f3$ax2 == 1]), "exclusion\n")
fam <- f3 %>% count(eje, familia) %>% group_by(eje) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>% arrange(eje, desc(n)) %>% ungroup()
print(as.data.frame(fam %>% group_by(eje) %>% slice_head(n = 8) %>% ungroup()))

cat("\nthe ladder — % of complete loci on the exclusion axis, by fraction:\n")
lad <- f3 %>% group_by(Frac) %>%
  summarise(ax1 = sum(ax2 == 0), ax2n = sum(ax2 == 1), n = n(),
            pct = 100 * mean(ax2), .groups = "drop") %>% arrange(match(Frac, FR_LEV))
print(as.data.frame(lad %>% mutate(pct = round(pct, 1))))
lad_life <- f3 %>% group_by(Frac, Life) %>%
  summarise(n = n(), pct = 100 * mean(ax2), .groups = "drop") %>%
  arrange(match(Frac, FR_LEV), Life)
print(as.data.frame(lad_life %>% mutate(pct = round(pct, 1))))

cat("\n  logistic, VIR as reference:\n")
m_frac <- glm(ax2 ~ Frac, data = f3, family = binomial)
s <- summary(m_frac)$coefficients
print(round(cbind(s[, c(1, 4)], OR = exp(s[, 1])), 4))
for (r in rownames(s)[-1]) {
  ci <- exp(s[r, 1] + c(-1.96, 1.96) * s[r, 2])
  add("axis composition", "Logistic regression", paste0(sub("^Frac", "", r), " vs VIR"),
      estimate = exp(s[r, 1]), ci_low = ci[1], ci_high = ci[2], p = s[r, 4])
}
cat("\n  the confound, tested both ways:\n")
cmh_report(with(f3, table(Life, ax2, Frac)), "lifestyle, stratified by fraction",
           "axis composition")
f3 <- f3 %>% mutate(Cellular = factor(ifelse(Frac == "VIR", "VIR", "cell-derived"),
                                      levels = c("VIR", "cell-derived")))
cmh_report(with(f3, table(Cellular, ax2, Life)),
           "cell-derived vs VIR, stratified by lifestyle", "axis composition")
m_life <- glm(ax2 ~ Life, data = f3, family = binomial)
m_both <- glm(ax2 ~ Frac + Life, data = f3, family = binomial)
m_full <- glm(ax2 ~ Frac + Life + loglen, data = f3, family = binomial)
for (z in list(list(m_frac, m_both, "lifestyle added to a model containing fraction"),
               list(m_life, m_both, "fraction added to a model containing lifestyle"),
               list(m_both, m_full, "log10(contig length) added to fraction + lifestyle"))) {
  a <- anova(z[[1]], z[[2]], test = "Chisq")
  cat(sprintf("  LR: %-52s p = %.4g\n", z[[3]], a$`Pr(>Chi)`[2]))
  add("axis composition", "Likelihood-ratio test", z[[3]], p = a$`Pr(>Chi)`[2])
}

# =============================================================================
rule("5. DEMOGRAPHY OF THE TWO AXES")
carr1 <- unique(f3$contig_id[f3$ax2 == 0]); carr2 <- unique(f3$contig_id[f3$ax2 == 1])
demo <- votu %>%
  mutate(Group = factor(case_when(
    vOTU_ID %in% carr1 & vOTU_ID %in% carr2 ~ "both axes",
    vOTU_ID %in% carr1 ~ "self-protection only",
    vOTU_ID %in% carr2 ~ "exclusion only",
    TRUE ~ "no complete locus"),
    levels = c("self-protection only", "exclusion only", "both axes", "no complete locus")))
print(as.data.frame(demo %>% group_by(Group) %>%
                      summarise(vOTUs = n(), CPM = round(sum(Total_TPM)),
                                pct_CPM = round(100 * sum(Total_TPM) / sum(demo$Total_TPM), 1),
                                median_CPM = round(median(Total_TPM), 1),
                                median_len_kb = round(median(contig_length) / 1000, 1),
                                .groups = "drop")))
w <- wilcox.test(Total_TPM ~ Group,
                 data = demo %>% filter(Group %in% c("self-protection only", "exclusion only")) %>%
                   droplevels())
cat(sprintf("\n  Wilcoxon, self-protection-only vs exclusion-only CPM: p = %.4g\n", w$p.value))
add("axis demography", "Wilcoxon rank-sum",
    "CPM per vOTU, self-protection-only vs exclusion-only", p = w$p.value)

# =============================================================================
rule("WRITE")
write_tsv(bind_rows(RES) %>% mutate(across(where(is.numeric), ~signif(.x, 4))),
          "stats_results_fig3.tsv")
write_tsv(BUT, "fig3_butterfly.tsv")
write_tsv(prev %>% left_join(cpm, by = c("Frac", "Life")), "fig3_panel_data.tsv")
write_tsv(lad %>% left_join(lad_life %>% rename(pct_life = pct, n_life = n), by = "Frac"),
          "fig3_axis_ladder.tsv")
write_tsv(fam, "fig3_families.tsv")
write_tsv(demo %>% select(vOTU_ID, Group, Total_TPM, contig_length), "fig3_demography.tsv")
cat("written. now run 05_Fig_3.R\n")
sink()
