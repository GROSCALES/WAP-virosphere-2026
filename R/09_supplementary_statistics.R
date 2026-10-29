# =============================================================================
# 09_supplementary_statistics.R
# Reproduces every test reported in Supplementary Tables S3 (Fig. 2) and S4
# (Fig. 3) from the three data files distributed with the paper:
#   Data S1 (vOTU catalogue), Data S2 (classified genes), Data S3 (loci),
#   either as distributed (.xlsx, sheet "Data") or exported to tab-separated text.
# Base R + MASS, plus readxl only when the inputs are .xlsx. Community-level
# tests on the eight libraries (Wilcoxon, PERMANOVA with exact p-values,
# betadisper) need the per-library abundance matrix and are run in 01_stats.R.
#
# WRITES two files:
#   supplementary_statistics.tsv             Table S3 and Table S4c
#   supplementary_statistics_categories.tsv  Tables S4a and S4b (four analysis
#                                            sets: main, >=10 kb or complete,
#                                            context score <=3, no proviruses)
# =============================================================================
suppressMessages(library(MASS))

args <- commandArgs(trailingOnly = TRUE)
guess <- function(stem) {                 # .xlsx as distributed, .tsv if exported
  cand <- paste0(stem, c(".xlsx", ".tsv", ".txt"))
  ok <- cand[file.exists(cand)]
  if (!length(ok)) stop("Not found: ", paste(cand, collapse = " / "), call. = FALSE)
  ok[1]
}
f_votu  <- if (length(args) >= 1) args[1] else guess("Data_S1_vOTU_catalogue")
f_genes <- if (length(args) >= 2) args[2] else guess("Data_S2_classified_genes")
f_loci  <- if (length(args) >= 3) args[3] else guess("Data_S3_defence_loci")
f_out   <- if (length(args) >= 4) args[4] else "supplementary_statistics.tsv"

# -----------------------------------------------------------------------------
# INPUT — Data S1-S3, as distributed (.xlsx, sheet "Data") or as plain text.
# .xlsx needs the readxl package; nothing else in this script does.
# -----------------------------------------------------------------------------
read_data <- function(path) {
  if (grepl("\\.xlsx?$", path, ignore.case = TRUE)) {
    if (!requireNamespace("readxl", quietly = TRUE))
      stop("Reading ", basename(path), " needs the readxl package.\n",
           "  install.packages(\"readxl\")  -- or save its 'Data' sheet as a\n",
           "  tab-separated .tsv and pass that instead.", call. = FALSE)
    d <- as.data.frame(readxl::read_excel(path, sheet = "Data"),
                       stringsAsFactors = FALSE)
    # read_excel returns NA for empty cells and read.delim returns ""; the
    # tests below compare strings, so the two are made identical here.
    d[] <- lapply(d, function(col) if (is.character(col))
      ifelse(is.na(col), "", col) else col)
    d
  } else {
    read.delim(path, stringsAsFactors = FALSE)
  }
}
v <- read_data(f_votu)
g <- read_data(f_genes)
l <- read_data(f_loci)
stopifnot(nrow(v) == 11100, nrow(g) == 3405, nrow(l) == 1718)

v$Frac   <- factor(v$Fraction, levels = c("VIR", "LIS", "BAC", "CUL"))
v$Life   <- factor(ifelse(grepl("^Temperate", v$Lifestyle), "Temperate", "Virulent"),
                   levels = c("Virulent", "Temperate"))
v$Orig   <- factor(ifelse(v$Origin == "Culture-dependent", "CD", "CI"),
                   levels = c("CI", "CD"))
v$len    <- v$Contig_length_bp
v$temp   <- v$Life == "Temperate"
v$Zone   <- factor(v$Station)

RES <- list()
add <- function(block, label, test, n = NA, est = NA, lo = NA, hi = NA, p = NA,
                q = NA, stat = NA, df = NA, note = "") {
  RES[[length(RES) + 1]] <<- data.frame(block, label, test, n, est, lo, hi, p, q,
                                        stat, df, note, stringsAsFactors = FALSE)
}
fish <- function(x, y) {           # x: factor (reference first); y: logical
  ft <- fisher.test(table(x, factor(y, levels = c(FALSE, TRUE))))
  c(est = unname(ft$estimate), lo = ft$conf.int[1], hi = ft$conf.int[2], p = ft$p.value)
}
mh <- function(x, y, s) {          # Mantel-Haenszel OR, RBG CI, CMH p (continuity-corrected)
  tb <- table(x, factor(y, levels = c(FALSE, TRUE)), droplevels(s))
  r  <- mantelhaen.test(tb, exact = FALSE)
  # Robins-Breslow-Greenland 95% CI, computed explicitly (mantelhaen.test
  # returns NA for some tables)
  a <- tb[1, 1, ]; b <- tb[1, 2, ]; c <- tb[2, 1, ]; d <- tb[2, 2, ]
  n <- a + b + c + d
  R <- a * d / n; S <- b * c / n; P <- (a + d) / n; Q <- (b + c) / n
  est <- sum(R) / sum(S)
  vr  <- sum(P * R) / (2 * sum(R)^2) + sum(P * S + Q * R) / (2 * sum(R) * sum(S)) +
         sum(Q * S) / (2 * sum(S)^2)
  ci  <- exp(log(est) + c(-1, 1) * qnorm(0.975) * sqrt(vr))
  if (!is.finite(est) || est == 0) ci <- c(NA, NA)
  stopifnot(isTRUE(all.equal(est, unname(r$estimate))) || !is.finite(est) || est == 0)
  c(est = est, lo = ci[1], hi = ci[2], p = r$p.value, stat = unname(r$statistic))
}
lr <- function(m0, m1) {           # likelihood-ratio test between nested models
  a <- anova(m0, m1, test = "Chisq")
  if ("LR stat." %in% names(a)) {  # negbin anova
    return(c(stat = a$`LR stat.`[2], df = a$`   df`[2], p = a$`Pr(Chi)`[2]))
  }
  c(stat = a$Deviance[2], df = a$Df[2], p = a$`Pr(>Chi)`[2])
}

# -----------------------------------------------------------------------------
# TABLE S3 — lifestyle (Fig. 2d and robustness)
# -----------------------------------------------------------------------------
r <- fish(v$Orig, v$temp)
add("S3_life", "Culture-dependent vs culture-independent", "Fisher",
    sprintf("%d vs %d", sum(v$Orig == "CD"), sum(v$Orig == "CI")), r["est"], r["lo"], r["hi"], r["p"])

cb <- v[v$Frac %in% c("BAC", "CUL"), ]; cb$Frac <- factor(cb$Frac, levels = c("BAC", "CUL"))
r <- fish(cb$Frac, cb$temp)
add("S3_life", "CUL vs BAC", "Fisher",
    sprintf("%d vs %d", sum(cb$Frac == "CUL"), sum(cb$Frac == "BAC")), r["est"], r["lo"], r["hi"], r["p"])

cb$lenbin <- cut(cb$len, c(0, 5e3, 1e4, 2e4, Inf))
r <- mh(cb$Frac, cb$temp, cb$lenbin)
add("S3_life", "CUL vs BAC, stratified by contig length", "CMH",
    "4 length strata", r["est"], r["lo"], r["hi"], r["p"], stat = r["stat"])

cmp <- v[!is.na(v$CheckV_completeness_pct), ]
cq  <- cmp[cmp$Frac %in% c("BAC", "CUL"), ]; cq$Frac <- factor(cq$Frac, levels = c("BAC", "CUL"))
cq$cbin <- cut(cq$CheckV_completeness_pct, quantile(cq$CheckV_completeness_pct, 0:4 / 4), include.lowest = TRUE)
r <- mh(cq$Frac, cq$temp, cq$cbin)
add("S3_life", "CUL vs BAC, stratified by CheckV completeness", "CMH",
    "4 completeness quartiles", r["est"], r["lo"], r["hi"], r["p"], stat = r["stat"],
    note = sprintf("%d vOTUs with a completeness estimate", nrow(cmp)))

hq <- v[v$CheckV_quality %in% c("Medium-quality", "High-quality", "Complete"), ]
r <- fish(hq$Orig, hq$temp)
add("S3_life", "Culture-dependent vs culture-independent, medium quality or better", "Fisher",
    sprintf("%d vs %d", sum(hq$Orig == "CD"), sum(hq$Orig == "CI")), r["est"], r["lo"], r["hi"], r["p"])

# Robustness: drop the 429 vOTUs with a SINGLE line of supporting evidence
# besides VirSorter2 (Supplementary Methods S5).
#
# CORRECTED. The earlier version of this block removed 329 vOTUs, split
# BAC 193 / VIR 128 / LIS 5 / CUL 3, and assumed all of them were virulent.
# Counted directly in votu_06_catalog.tsv (Passes_catalog == "Yes",
# N_evidence == 1) the set is 429, split BAC 258 / VIR 157 / LIS 10 / CUL 4,
# and 29 of them are TEMPERATE (BAC 19, VIR 8, LIS 2), so both the temperate
# and the virulent counts have to be reduced. The old block reduced only the
# virulent ones.
#
# ON THE NAME. VirSorter2 is the baseline: every catalogue sequence was called
# by it, and N_evidence counts the corroborating tools, which is why
# Evidence_sources never mentions VirSorter2. "Supported by VirSorter2 alone"
# would be N_evidence == 0, and there are none: all 11,100 carry
# Has_evidence == "Yes". The 429 removed here are those with exactly one
# corroborating tool (419 PhaBOX, 10 Provirus).
#
# The counts are written out because Data S1, which is all this script reads,
# does not carry the evidence column. One line reproduces them from
# votu_06_catalog.tsv:
#   s <- t[t$Passes_catalog == "Yes" & t$N_evidence == 1, ]
#   table(s$Fraction, grepl("^Temperate", s$Lifestyle))
VS2_T <- c(BAC =  19, VIR =   8, LIS =  2, CUL = 0)   # temperate
VS2_V <- c(BAC = 239, VIR = 149, LIS =  8, CUL = 4)   # virulent
stopifnot(sum(VS2_T) + sum(VS2_V) == 429)

tt <- table(v$Frac, v$temp)
cnt <- function(fr, t) tt[fr, as.character(t)]
dep_t <- cnt("CUL", TRUE)  + cnt("LIS", TRUE)  - VS2_T["CUL"] - VS2_T["LIS"]
dep_v <- cnt("CUL", FALSE) + cnt("LIS", FALSE) - VS2_V["CUL"] - VS2_V["LIS"]
ind_t <- cnt("BAC", TRUE)  + cnt("VIR", TRUE)  - VS2_T["BAC"] - VS2_T["VIR"]
ind_v <- cnt("BAC", FALSE) + cnt("VIR", FALSE) - VS2_V["BAC"] - VS2_V["VIR"]
stopifnot(dep_t + dep_v + ind_t + ind_v == 11100 - 429)
ft <- fisher.test(matrix(c(ind_v, ind_t, dep_v, dep_t), 2, byrow = TRUE))
add("S3_life", "Culture-dependent vs culture-independent, without single-evidence vOTUs", "Fisher",
    sprintf("%d vs %d", dep_t + dep_v, ind_t + ind_v), ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value)
ft <- fisher.test(matrix(c(cnt("BAC", FALSE) - VS2_V["BAC"], cnt("BAC", TRUE) - VS2_T["BAC"],
                           cnt("CUL", FALSE) - VS2_V["CUL"], cnt("CUL", TRUE) - VS2_T["CUL"]), 2, byrow = TRUE))
add("S3_life", "CUL vs BAC, without single-evidence vOTUs", "Fisher",
    sprintf("%d vs %d", sum(tt["CUL", ]) - VS2_T["CUL"] - VS2_V["CUL"],
                        sum(tt["BAC", ]) - VS2_T["BAC"] - VS2_V["BAC"]),
    ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value)

vb <- v[v$Frac %in% c("BAC", "VIR"), ]; vb$Frac <- factor(vb$Frac, levels = c("BAC", "VIR"))
r <- fish(vb$Frac, vb$temp)
add("S3_proc", "VIR vs BAC", "Fisher",
    sprintf("%d vs %d", sum(vb$Frac == "VIR"), sum(vb$Frac == "BAC")), r["est"], r["lo"], r["hi"], r["p"])
lc <- v[v$Frac %in% c("CUL", "LIS"), ]; lc$Frac <- factor(lc$Frac, levels = c("CUL", "LIS"))
r <- fish(lc$Frac, lc$temp)
add("S3_proc", "LIS vs CUL", "Fisher",
    sprintf("%d vs %d", sum(lc$Frac == "LIS"), sum(lc$Frac == "CUL")), r["est"], r["lo"], r["hi"], r["p"])

# -----------------------------------------------------------------------------
# TABLE S4a — category-level carriage (Fig. 3a-b)
# -----------------------------------------------------------------------------
fig3 <- l[l$In_Figure3 == "yes", ]
stopifnot(nrow(fig3) == 773)

# context score of a locus: every defence gene inside its coordinates <= 3
gd <- g[g$Gene_class == "Defence", c("vOTU_ID", "Start_bp", "End_bp", "Context_score")]
names(gd)[2:3] <- c("gs", "ge")
mm <- merge(fig3[, c("Locus_ID", "vOTU_ID", "Start_bp", "End_bp")], gd, by = "vOTU_ID")
mm <- mm[mm$gs >= mm$Start_bp & mm$ge <= mm$End_bp, ]
le3 <- tapply(mm$Context_score <= 3, mm$Locus_ID, all)
fig3$le3 <- unname(le3[fig3$Locus_ID])
stopifnot(!any(is.na(fig3$le3)))

amg <- g[g$Gene_class == "AMG", ]
DEF_CATS <- c("Genome modification", "Nucleic acid degradation",
              "Abortive infection", "Other / uncharacterised systems")
AMG_CATS <- names(sort(tapply(amg$vOTU_ID, amg$Category, function(x) length(unique(x))),
                       decreasing = TRUE))

carriers <- function(klass, cat, scen) {
  if (klass == "Defence") {
    d <- fig3[fig3$Category == cat, ]
    if (scen == "le3") d <- d[d$le3, ]
    unique(d$vOTU_ID)
  } else {
    d <- amg[amg$Category == cat, ]
    if (scen == "le3") d <- d[d$Context_score <= 3, ]
    unique(d$vOTU_ID)
  }
}
SCEN <- list(main  = v,
             kb10  = v[v$Subset_10kb_or_complete == "Yes", ],
             le3   = v,
             nopro = v[v$Lifestyle_source != "CheckV provirus", ])
stopifnot(nrow(SCEN$kb10) == 3810, nrow(SCEN$nopro) == 10729)

CAT <- list()
for (klass in c("Defence", "AMG")) {
  cats <- if (klass == "Defence") DEF_CATS else AMG_CATS
  for (sc in names(SCEN)) {
    d <- SCEN[[sc]]
    for (ct in cats) {
      has <- d$vOTU_ID %in% carriers(klass, ct, sc)
      o <- if (sc == "nopro") rep(NA, 5) else mh(d$Orig, has, d$Life)
      f <- mh(d$Life, has, d$Frac)
      unit_n <- if (klass == "Defence") sum(fig3$Category == ct) else sum(amg$Category == ct)
      CAT[[length(CAT) + 1]] <- data.frame(
        class = klass, scen = sc, category = ct, n_units = unit_n,
        carriers = sum(has), cd = sum(has & d$Orig == "CD"), ci = sum(has & d$Orig == "CI"),
        tp = sum(has & d$Life == "Temperate"), vr = sum(has & d$Life == "Virulent"),
        or_o = o[1], lo_o = o[2], hi_o = o[3], p_o = o[4],
        or_l = f[1], lo_l = f[2], hi_l = f[3], p_l = f[4], stringsAsFactors = FALSE)
    }
  }
}
CAT <- do.call(rbind, CAT)
# BH within gene class x contrast x analysis set
CAT$q_o <- NA; CAT$q_l <- NA
for (k in unique(CAT$class)) for (s in unique(CAT$scen)) {
  i <- CAT$class == k & CAT$scen == s
  if (s != "nopro") CAT$q_o[i] <- p.adjust(CAT$p_o[i], "BH")
  CAT$q_l[i] <- p.adjust(CAT$p_l[i], "BH")
}
write.table(CAT, sub("\\.tsv$", "_categories.tsv", f_out), sep = "\t", quote = FALSE, row.names = FALSE)

# -----------------------------------------------------------------------------
# TABLE S4b — other tests behind Fig. 3
# -----------------------------------------------------------------------------
v$hasD <- v$N_defence_genes > 0
v$hasA <- v$N_AMGs > 0
r <- mh(v$Life, v$hasD, v$Frac)
add("S4_carry", "Defence gene ~ lifestyle | fraction", "MH", nrow(v), r["est"], r["lo"], r["hi"], r["p"], stat = r["stat"])
r <- mh(v$Orig, v$hasD, v$Life)
add("S4_carry", "Defence gene ~ origin | lifestyle", "MH", nrow(v), r["est"], r["lo"], r["hi"], r["p"], stat = r["stat"])
r <- mh(v$Life, v$hasA, v$Frac)
add("S4_carry", "AMG ~ lifestyle | fraction", "MH", nrow(v), r["est"], r["lo"], r["hi"], r["p"], stat = r["stat"])
r <- mh(v$Orig, v$hasA, v$Life)
add("S4_carry", "AMG ~ origin | lifestyle", "MH", nrow(v), r["est"], r["lo"], r["hi"], r["p"], stat = r["stat"])
np <- SCEN$nopro; np$hasD <- np$N_defence_genes > 0
r <- mh(np$Life, np$hasD, np$Frac)
add("S4_carry", "Defence gene ~ lifestyle | fraction, without CheckV proviruses", "MH",
    nrow(np), r["est"], r["lo"], r["hi"], r["p"], stat = r["stat"])

nl <- table(fig3$vOTU_ID)
v$n_loci <- as.integer(ifelse(v$vOTU_ID %in% names(nl), nl[v$vOTU_ID], 0))
stopifnot(sum(v$n_loci) == 773)
for (y in c("n_loci", "N_AMGs")) {
  full <- glm.nb(as.formula(paste(y, "~ Orig * Life + Zone + log(len)")), data = v)
  noO  <- glm.nb(as.formula(paste(y, "~ Life + Zone + log(len)")), data = v)
  noL  <- glm.nb(as.formula(paste(y, "~ Orig + Zone + log(len)")), data = v)
  a <- lr(noO, full); add("S4_nb", paste(y, ": origin, given lifestyle"), "NB LR", nrow(v), p = a["p"], stat = a["stat"], df = a["df"])
  a <- lr(noL, full); add("S4_nb", paste(y, ": lifestyle, given origin"), "NB LR", nrow(v), p = a["p"], stat = a["stat"], df = a["df"])
}

f3 <- merge(fig3, v[, c("vOTU_ID", "Frac", "Life", "len")], by = "vOTU_ID")
f3$ax2 <- as.integer(f3$Axis == "Exclusion")
stopifnot(nrow(f3) == 773, sum(f3$ax2) == 174)
m_frac <- glm(ax2 ~ Frac, data = f3, family = binomial)
s <- summary(m_frac)$coefficients
for (rn in rownames(s)[-1]) {
  add("S4_axis", paste(sub("^Frac", "", rn), "vs VIR"), "Logistic",
      sum(f3$Frac %in% c("VIR", sub("^Frac", "", rn))),
      exp(s[rn, 1]), exp(s[rn, 1] - 1.96 * s[rn, 2]), exp(s[rn, 1] + 1.96 * s[rn, 2]), s[rn, 4])
}
r <- mh(f3$Life, f3$ax2 == 1, f3$Frac)
add("S4_axis", "Lifestyle, within fractions", "CMH", nrow(f3), r["est"], r["lo"], r["hi"], r["p"], stat = r["stat"])
f3$Cell <- factor(ifelse(f3$Frac == "VIR", "VIR", "cell"), levels = c("VIR", "cell"))
r <- mh(f3$Cell, f3$ax2 == 1, f3$Life)
add("S4_axis", "Cell-derived vs VIR, within lifestyles", "CMH", nrow(f3), r["est"], r["lo"], r["hi"], r["p"], stat = r["stat"])
m_life <- glm(ax2 ~ Life, data = f3, family = binomial)
m_both <- glm(ax2 ~ Frac + Life, data = f3, family = binomial)
m_full <- glm(ax2 ~ Frac + Life + log10(len), data = f3, family = binomial)
a <- lr(m_frac, m_both); add("S4_axis", "Lifestyle added to fraction", "LR", nrow(f3), p = a["p"], stat = a["stat"], df = a["df"])
a <- lr(m_life, m_both); add("S4_axis", "Fraction added to lifestyle", "LR", nrow(f3), p = a["p"], stat = a["stat"], df = a["df"])
a <- lr(m_both, m_full); add("S4_axis", "Contig length added to both", "LR", nrow(f3), p = a["p"], stat = a["stat"], df = a["df"])

ft <- fisher.test(table(v$hasA, v$hasD))
add("S4_cooc", "AMG and defence gene in the same vOTU, crude", "Fisher", nrow(v),
    ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value)
mco <- glm(hasD ~ hasA + log(len) + Orig, data = v, family = binomial)
s <- summary(mco)$coefficients; ci <- suppressMessages(confint(mco))
add("S4_cooc", "AMG and defence gene in the same vOTU, adjusted for length and origin", "Logistic",
    nrow(v), exp(s["hasATRUE", 1]), exp(ci["hasATRUE", 1]), exp(ci["hasATRUE", 2]), s["hasATRUE", 4])

out <- do.call(rbind, RES)
write.table(out, f_out, sep = "\t", quote = FALSE, row.names = FALSE)
print(out[, c("block", "label", "n", "est", "lo", "hi", "p", "stat", "df")], digits = 4)
