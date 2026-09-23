# =============================================================================
# palette_maps.R — one palette and one vocabulary for the genome-map figures
# (Figures 4, 5 and 6). Sourced by all three; no map script defines its own.
#
# FOUR THINGS CHANGED, ALL OF THEM VOCABULARY RATHER THAN COLOUR
#
#   1. "Defense" -> "Defence". The maps were the only figures in the paper
#      using the American spelling, and it appeared six times in a legend
#      printed next to Figure 3, which spells it the other way.
#
#   2. "counter-defence" -> "Self-protection". The two axes are now called
#      Self-protection and Exclusion (see palette_paper.R): "Defence" is the
#      name of the gene CLASS, so using it for one of the two axes made the
#      same word mean two things at two levels.
#
#   3. "AMG: sulfur relay (...)" -> "AMG: Folding, sorting and degradation".
#      Its 20 genes really are sulfur-relay enzymes (16 mec, 2 moeB, 2 tusA)
#      and KEGG files them under protein folding, which is misleading -- but
#      renaming a published category inside a figure legend asserts by decree
#      something that belongs in the Discussion, where it can be argued. Every
#      figure and Table S3 now use the same 19 category names.
#
#   4. The error message in relabel_map() was in Spanish.
#
# Colours are untouched.
# =============================================================================

MAP_CATEGORIES <- c(
  # --- context -------------------------------------------------------------
  "Host (MAG / flanking)",
  "Structural / viral assembly",
  "Hypothetical / unannotated",
  
  # --- lysogeny lifecycle (NOT defence; axes 3A/3B in this study) ----------
  "Integration machinery",
  "SOS response / prophage induction",
  
  # --- regulators that are NOT defence -------------------------------------
  "Transcriptional regulator",
  
  # --- defence, prefixes matching Table S3 --------------------------------
  "Defence: Nucleic acid degradation",
  "Defence: Genome modification",
  "Defence: Toxin-antitoxin (Abortive infection)",
  "Defence: Anti-restriction (self-protection)",
  "Defence: Other / uncharacterised",
  
  # --- diversity-generating ------------------------------------------------
  "Hypermutation (DGR)",
  
  # --- auxiliary metabolic genes, prefixes matching Table S3 --------------
  # Energy metabolism was missing. Figure 6 needs it for cysH (PAPS reductase,
  # sulfur assimilation -> KEGG Energy metabolism), which was being drawn in
  # the cofactor/vitamin colour purely because there was no entry for it. That
  # put a category name in a figure legend that contradicts Table S3, which is
  # the one thing a shared palette exists to prevent.
  "AMG: Energy metabolism",
  "AMG: Metabolism of cofactors and vitamins",
  "AMG: Carbohydrate metabolism",
  "AMG: Folding, sorting and degradation"
)

MAP_PALETTE <- c(
  "Host (MAG / flanking)"                          = "#4D4D4D",
  "Structural / viral assembly"                    = "#377EB8",
  "Hypothetical / unannotated"                     = "#F0F0F0",
  
  # lysogeny lifecycle
  "Integration machinery"                          = "#E7298A",
  "SOS response / prophage induction"              = "#F781BF",
  
  # defence
  "Defence: Nucleic acid degradation"              = "#984EA3",
  "Defence: Genome modification"                   = "#E41A1C",
  "Defence: Toxin-antitoxin (Abortive infection)"  = "#FF7F00",
  "Transcriptional regulator"                      = "#8C564B",
  "Defence: Anti-restriction (self-protection)"    = "#005F73",
  # Added for Figure 6. "Other / uncharacterised" is 386 genes and 56 complete
  # loci in the catalogue, so the maps needed it. A mid grey on purpose: it
  # sits between the dark host grey and the near-white hypothetical, and the
  # neutrality is the message -- these are systems whose mechanism is not
  # known, not a mechanism of their own.
  "Defence: Other / uncharacterised"               = "#9E9E9E",
  
  "Hypermutation (DGR)"                            = "#00BFFF",
  
  # AMGs
  "AMG: Energy metabolism"                         = "#66C2A5",
  "AMG: Metabolism of cofactors and vitamins"      = "#FDBF6F",
  "AMG: Carbohydrate metabolism"                   = "#B2DF8A",
  "AMG: Folding, sorting and degradation"          = "#E3B62C"
)

stopifnot(setequal(MAP_CATEGORIES, names(MAP_PALETTE)))
stopifnot(!any(duplicated(MAP_PALETTE)))

# Old spellings are kept as inputs so the map scripts do not have to change:
# whatever they pass in, it comes out in the current vocabulary.
MAP_RELABEL <- c(
  # --- the mechanism-first labels tried in Figure 5 ------------------------
  # "Restriction-modification (Nucleic acid degradation)" reads better than
  # "Defence: Nucleic acid degradation", but it cannot live in a shared
  # palette: naming the mechanism commits every gene of that colour to it,
  # and Figure 4's purple gene is `orn`, an oligoribonuclease, not a
  # restriction enzyme. Figure 5's purple gene is mcrA, where it does fit.
  # Category names are correct for any gene in the category, so those are
  # what the shared file carries; the mechanism goes in the figure legend
  # text, where it can be specific to the gene actually shown.
  "Restriction-modification (Nucleic acid degradation)" = "Defence: Nucleic acid degradation",
  "DNA methylation (Genome modification)"   = "Defence: Genome modification",
  "Toxin-antitoxin (Abortive infection)"    = "Defence: Toxin-antitoxin (Abortive infection)",
  "Anti-restriction (Self-protection)"      = "Defence: Anti-restriction (self-protection)",
  "AVG: sulfur relay (Folding, sorting & degradation)" = "AMG: Folding, sorting and degradation",
  "AVG: Folding, sorting and degradation"   = "AMG: Folding, sorting and degradation",
  
  # --- current names, and the American spellings they replace -------------
  "Defense: Restriction Modification"       = "Defence: Nucleic acid degradation",
  "Defense: Nucleic acid degradation"       = "Defence: Nucleic acid degradation",
  "Defense: DNA Methylation"                = "Defence: Genome modification",
  "Defense: Genome modification"            = "Defence: Genome modification",
  "Defense: Toxin-Antitoxin"                = "Defence: Toxin-antitoxin (Abortive infection)",
  "Defense: Toxin-antitoxin (Abortive infection)" = "Defence: Toxin-antitoxin (Abortive infection)",
  "Defense: Phage Communication"            = "Transcriptional regulator",
  "Defense: Anti-restriction"               = "Defence: Anti-restriction (self-protection)",
  "Defense: Anti-restriction (counter-defence)"   = "Defence: Anti-restriction (self-protection)",
  "Anti-restriction (counter-defence)"      = "Defence: Anti-restriction (self-protection)",
  "Defense: Innate Immune System"           = "SOS response / prophage induction",
  
  # --- capitalisation variants --------------------------------------------
  "Host (MAG / Flanking)"                   = "Host (MAG / flanking)",
  "Structural / Viral Assembly"             = "Structural / viral assembly",
  "Hypothetical / Unannotated"              = "Hypothetical / unannotated",
  "Integration Machinery"                   = "Integration machinery",
  "Hypermutation (DGR)"                     = "Hypermutation (DGR)",
  
  # --- AMG shorthands and the sulfur-relay label this file used to invent --
  "AMG: Cofactors & Vitamins"               = "AMG: Metabolism of cofactors and vitamins",
  "AMG: Carbohydrate & Glycan"              = "AMG: Carbohydrate metabolism",
  "AMG: sulfur relay (Folding, sorting & degradation)" = "AMG: Folding, sorting and degradation",
  "AMG: Folding, sorting and degradation"   = "AMG: Folding, sorting and degradation",
  
  # --- bare category names, as they appear in Table S3 --------------------
  "Nucleic acid degradation"                = "Defence: Nucleic acid degradation",
  "Genome modification"                     = "Defence: Genome modification",
  "Other / uncharacterised systems"         = "Defence: Other / uncharacterised",
  "Abortive infection"                      = "Defence: Toxin-antitoxin (Abortive infection)",
  "Energy metabolism"                       = "AMG: Energy metabolism",
  "Metabolism of cofactors and vitamins"    = "AMG: Metabolism of cofactors and vitamins",
  "Carbohydrate metabolism"                 = "AMG: Carbohydrate metabolism",
  "Folding, sorting and degradation"        = "AMG: Folding, sorting and degradation"
)

relabel_map <- function(x) {
  out <- ifelse(x %in% names(MAP_RELABEL), unname(MAP_RELABEL[x]), x)
  unknown <- setdiff(unique(out), MAP_CATEGORIES)
  if (length(unknown)) {
    stop("Category not covered by palette_maps.R: ",
         paste(unknown, collapse = " | "),
         "\n  Add it to MAP_CATEGORIES and MAP_PALETTE, or map it in ",
         "MAP_RELABEL.", call. = FALSE)
  }
  present <- MAP_CATEGORIES[MAP_CATEGORIES %in% out]
  factor(out, levels = present)
}
