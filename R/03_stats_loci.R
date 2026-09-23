# =============================================================================
# 03_stats_loci.R — the Defence side of Figure 3, counted on LOCI.
#
# VERSION 23/9. THREE CHANGES over the 28/8 script, everything else identical:
#
#   1. THE CATEGORY ENRICHMENT IS NOW MANTEL-HAENSZEL, NOT CRUDE FISHER.
#      Origin and lifestyle are confounded in this design (CUL loci are 91%
#      temperate, VIR loci 27%), so a crude odds ratio for one of them carries
#      part of the other. Each contrast is now estimated with the other held
#      constant, which is what Table S4a reports and what Figure 3b draws:
#        origin    = culture-dependent vs -independent, stratified by lifestyle
#        lifestyle = temperate vs virulent, stratified by fraction (4 strata)
#      The crude Fisher odds ratios are still written, under the section name
#      "category enrichment (loci, unadjusted)", so the two can be compared;
#      the figure ignores that section.
#
#   2. NO ABSOLUTE PATH. Run with the working directory set to the folder that
#      holds the tables. Input names are tried with and without the "_28_8"
#      suffix and the file actually opened is printed.
#
#   3. Integer-overflow guard in mantelhaen.test: the 2x2xK tables reach
#      11,100 vOTUs and R's CMH does its arithmetic in integers, which
#      silently returns NA for the confidence interval. The tables are cast to
#      double before the test.
#
# RUN ORDER   01_stats.R -> 03_stats_loci.R -> 04_stats_fig3.R -> 05_Fig_3.R
# OUTPUTS     07_long_table_loci_for_fig3.tsv, votu_n_defence_loci.tsv,
#             stats_results_loci_defence.tsv
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(MASS)
})
select <- dplyr::select   # MASS::select shadows dplyr::select; see 01_stats.R

pick <- function(...) {
  f <- c(...); ok <- f[file.exists(f)]
  if (!length(ok)) stop("No file found: ", paste(f, collapse = " / "))
  cat("   input:", ok[1], "\n"); ok[1]
}
F_LOCI   <- pick("10_loci_defensa_28_8.tsv", "10_loci_defensa.tsv")
F_VOTU   <- pick("07_votu_with_counts_28_8.tsv", "07_votu_with_counts.tsv")
F_CATMAP <- pick("06_category_map_28_8.tsv", "06_category_map.tsv")

results <- tibble(section = character(), test = character(),
                  comparison = character(), n1 = double(), n2 = double(),
                  estimate = double(), ci_low = double(), ci_high = double(),
                  p = double(), q = double(), note = character())
add <- function(section, test, comparison, n1 = NA, n2 = NA, estimate = NA,
                ci_low = NA, ci_high = NA, p = NA, q = NA, note = "") {
  results <<- bind_rows(results, tibble(
    section = section, test = test, comparison = comparison,
    n1 = n1, n2 = n2, estimate = estimate, ci_low = ci_low, ci_high = ci_high,
    p = p, q = q, note = note))
}

# =============================================================================
# READ
# =============================================================================
loci <- read.delim(F_LOCI, stringsAsFactors = FALSE)
votu <- read.delim(F_VOTU, stringsAsFactors = FALSE) %>%
  filter(Passes_catalog == "Yes") %>%
  mutate(
    Life_Gen = factor(ifelse(grepl("^Temperate", Lifestyle), "Temperate", "Virulent"),
                      levels = c("Virulent", "Temperate")),
    Fraction_group = factor(Fraction_group,
                            levels = c("Culture-independent", "Culture-dependent")),
    # the four strata of the lifestyle contrast, ordered from free particles
    # to cultured cells, as in Figure 3e
    Frac = factor(c(VIRFC = "VIR", LIS = "LIS", CB22 = "BAC", GLI = "CUL")[Fraction],
                  levels = c("VIR", "LIS", "BAC", "CUL"))
  )
stopifnot(!any(is.na(votu$Frac)))

def_loci <- loci %>%
  filter(completitud_locus == "completo",
         eje %in% c("1 Contra-defensa", "2 Defensa"))
cat("Completed defence loci kept:", nrow(def_loci), "on",
    n_distinct(def_loci$contig_id), "vOTUs\n")
cat("(expected after the distinct-gene rule: 599 + 174 = 773)\n")

# ---------------------------------------------------------------------------
# CATEGORY LABELS: the published English six, NOT mapa_v5's `subcategoria`.
# Unchanged from the 28/8 script.
# ---------------------------------------------------------------------------
catmap <- read.delim(F_CATMAP, stringsAsFactors = FALSE)
stopifnot(all(c("category_raw", "category") %in% names(catmap)))
catmap <- catmap %>% select(category_raw, category) %>% distinct()
cat("Category map:", nrow(catmap), "raw labels ->",
    n_distinct(catmap$category), "published categories\n")
lookup <- setNames(catmap$category, catmap$category_raw)

CAT_PRIORITY <- c(
  "CRISPR-Cas",
  "Abortive infection",
  "Nucleic acid degradation",
  "Genome modification",
  "Host manipulation (phage-encoded)",
  "Other / uncharacterised systems"
)
ambiguous <- character(0)
resolve_cat <- function(g) {
  parts <- trimws(strsplit(g, ";")[[1]])
  cats  <- unique(unname(lookup[parts]))
  cats  <- cats[!is.na(cats)]
  if (length(cats) == 0) return(NA_character_)
  if (length(cats) > 1) {
    ambiguous <<- c(ambiguous, g)
    return(cats[order(match(cats, CAT_PRIORITY))][1])
  }
  cats[1]
}
def_loci <- def_loci %>% mutate(category_en = vapply(genes, resolve_cat, character(1)))
if (length(ambiguous) > 0)
  cat("\n", length(ambiguous), " locus/loci span two categories; resolved by CAT_PRIORITY\n", sep = "")

# The category must follow the reassignment of rescinded RM systems ...
REASSIGN_CATEGORY <- c(
  "RM incompleto -> MTasa huerfana"    = "Genome modification",
  "RM incompleto -> nucleasa putativa" = "Nucleic acid degradation"
)
if ("reasignacion" %in% names(def_loci)) {
  def_loci <- def_loci %>%
    mutate(category_en = ifelse(
      !is.na(reasignacion) & reasignacion %in% names(REASSIGN_CATEGORY),
      unname(REASSIGN_CATEGORY[reasignacion]), category_en))
} else {
  stop("no `reasignacion` column in ", F_LOCI)
}
# ... and the same rule tied to the family, for the orphan methyltransferases
# that never went through an RM rescission.
if ("familia" %in% names(def_loci)) {
  def_loci <- def_loci %>%
    mutate(category_en = ifelse(familia == "MTasa huerfana",
                                "Genome modification", category_en))
}
cat("\nCategory x axis, after both overrides:\n")
print(as.data.frame(def_loci %>%
                      mutate(Axis = ifelse(eje == "1 Contra-defensa", "Self-protection", "Exclusion")) %>%
                      count(category_en, Axis) %>%
                      pivot_wider(names_from = Axis, values_from = n, values_fill = 0)))
stopifnot(!any(is.na(def_loci$category_en)))

# =============================================================================
# 1. n_Defence_loci per vOTU
# =============================================================================
n_per_votu <- def_loci %>% count(contig_id, name = "n_Defence_loci")
votu <- votu %>%
  left_join(n_per_votu, by = c("vOTU_ID" = "contig_id")) %>%
  mutate(n_Defence_loci = replace_na(n_Defence_loci, 0))
write.table(votu %>% select(vOTU_ID, n_Defence_loci),
            "votu_n_defence_loci.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# =============================================================================
# 2. long format, ONE ROW PER COMPLETED LOCUS
# =============================================================================
tpm_cols <- grep("^TPM_(DE|EL|MB|OR)_(CB22|GLI)$", names(votu), value = TRUE)
long_loci <- def_loci %>%
  transmute(vOTU_ID = contig_id, category_final = category_en,
            gene_class = "Defence") %>%
  left_join(votu %>% select(vOTU_ID, Fraction_group, Life_Gen, contig_length,
                            all_of(tpm_cols)), by = "vOTU_ID")
write.table(long_loci, "07_long_table_loci_for_fig3.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("\nCompleted loci by published category (axes 1+2 only):\n")
print(as.data.frame(def_loci %>% count(category_en, sort = TRUE)))

# =============================================================================
# 3. PER-CATEGORY ENRICHMENT — Mantel-Haenszel, one contrast per split
#    Carriage = the vOTU has >=1 completed locus of that category.
#    BH within each split, because the two answer different questions.
# =============================================================================
carriers <- long_loci %>% distinct(vOTU_ID, category_final)

cmh <- function(has, exposure, strata, lv_lo, lv_hi) {
  tb <- table(factor(exposure, levels = c(lv_lo, lv_hi)), has, strata)
  tb <- array(as.double(tb), dim(tb), dimnames(tb))   # see header, change 3
  r  <- mantelhaen.test(tb, exact = FALSE)
  c(OR = unname(r$estimate), lo = r$conf.int[1], hi = r$conf.int[2], p = r$p.value)
}

per_cat_loci <- function(var, strat, lv_lo, lv_hi, tag, strat_label) {
  res <- map_dfr(sort(unique(carriers$category_final)), function(cc) {
    d <- votu %>% mutate(has = vOTU_ID %in% carriers$vOTU_ID[carriers$category_final == cc])
    mh <- cmh(d$has, d[[var]], d[[strat]], lv_lo, lv_hi)
    fi <- fisher.test(table(factor(d[[var]], levels = c(lv_lo, lv_hi)), d$has))
    tibble(category = cc,
           n_hi = sum(d$has & d[[var]] == lv_hi),
           n_lo = sum(d$has & d[[var]] == lv_lo),
           OR = mh["OR"], lo = mh["lo"], hi = mh["hi"], p = mh["p"],
           OR_crude = unname(fi$estimate), p_crude = fi$p.value)
  }) %>%
    mutate(q = p.adjust(p, method = "BH"),
           q_crude = p.adjust(p_crude, method = "BH")) %>%
    arrange(desc(OR))
  cat("\nDefence (complete loci), ", lv_hi, " vs ", lv_lo,
      ", stratified by ", strat_label, ":\n", sep = "")
  print(as.data.frame(res %>% mutate(across(where(is.numeric), ~signif(.x, 4)))))
  for (i in seq_len(nrow(res))) {
    add("category enrichment (loci)", "Cochran-Mantel-Haenszel + BH",
        paste0(res$category[i], " ~ ", tag),
        n1 = res$n_hi[i], n2 = res$n_lo[i], estimate = res$OR[i],
        ci_low = res$lo[i], ci_high = res$hi[i], p = res$p[i], q = res$q[i],
        note = paste0("carriage of >=1 complete locus; stratified by ",
                      strat_label, "; BH q = ", signif(res$q[i], 4)))
    add("category enrichment (loci, unadjusted)", "Fisher exact + BH",
        paste0(res$category[i], " ~ ", tag),
        n1 = res$n_hi[i], n2 = res$n_lo[i], estimate = res$OR_crude[i],
        p = res$p_crude[i], q = res$q_crude[i],
        note = "crude; kept for comparison, not used in Fig. 3")
  }
  invisible(res)
}
per_cat_loci("Fraction_group", "Life_Gen", "Culture-independent", "Culture-dependent",
             "Origin", "lifestyle")
per_cat_loci("Life_Gen", "Frac", "Virulent", "Temperate",
             "Life", "fraction (4 strata)")

# =============================================================================
# 4. fraction x lifestyle on the NUMBER of complete loci per vOTU
# =============================================================================
m <- glm.nb(n_Defence_loci ~ Fraction_group * Life_Gen + Zone + log(contig_length),
            data = votu)
cat("\nn_Defence_loci ~ Fraction_group * Life_Gen + Zone + log(length)\n")
m0 <- glm.nb(n_Defence_loci ~ Life_Gen + Zone + log(contig_length), data = votu)
an <- anova(m0, m)
cat("  adding origin to a model that already has lifestyle:  LR p =",
    signif(an$`Pr(Chi)`[2], 4), "\n")
add("fraction vs lifestyle (loci)", "Likelihood-ratio test",
    "n_Defence_loci: fraction beyond lifestyle", p = an$`Pr(Chi)`[2],
    note = "NB GLM, adjusted for length and station; complete loci, axes 1-2")
m1 <- glm.nb(n_Defence_loci ~ Fraction_group + Zone + log(contig_length), data = votu)
an1 <- anova(m1, m)
cat("  adding lifestyle to a model that already has origin:  LR p =",
    signif(an1$`Pr(Chi)`[2], 4), "\n")
add("fraction vs lifestyle (loci)", "Likelihood-ratio test",
    "n_Defence_loci: lifestyle beyond fraction", p = an1$`Pr(Chi)`[2],
    note = "NB GLM, adjusted for length and station; complete loci, axes 1-2")

results <- results %>%
  mutate(sig_basis = coalesce(q, p),
         signif = case_when(is.na(sig_basis) ~ "", sig_basis < 0.001 ~ "***",
                            sig_basis < 0.01 ~ "**", sig_basis < 0.05 ~ "*",
                            TRUE ~ "ns")) %>%
  select(-sig_basis)

write.table(results, "stats_results_loci_defence.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("\nWrote: 07_long_table_loci_for_fig3.tsv, votu_n_defence_loci.tsv,",
    "stats_results_loci_defence.tsv\n")
cat("Now run 04_stats_fig3.R, then 05_Fig_3.R\n")
