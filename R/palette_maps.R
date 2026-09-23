# =============================================================================
# palette_maps.R — one palette and one vocabulary for the genome-map figures
# (Figures 4, 5 and 6). Sourced by all three; no map script defines its own.
#
# VERSION 23/9 (second pass). What changed and why:
#
#   1. NEW COLOURS. The old set mixed three warm tones that meant different
#      things in different figures (AMG cofactors, AMG folding, toxin-antitoxin)
#      and used fully saturated primaries. The new set is built from Paul Tol's
#      qualitative schemes, which are designed for colour-vision deficiency and
#      for print, and it is organised by family so that a reader learns it once:
#
#        neutral        context: host, other phage genes, unannotated
#        red / wine     defence
#        yellow / green defence-unrelated metabolism (AMG)
#        magenta, cyan  the two lifecycle accents (integration, DGR)
#
#      Every within-figure pair was checked for colour-blind separation. The
#      one pair below the usual threshold is "Defence: Other / uncharacterised"
#      (grey) against "Defence: Genome modification" (red) for deuteranopes;
#      both categories are always directly labelled in the figures, which is
#      the accepted relief for that case.
#
#   2. "Hypermutation (DGR)" is now "Diversity-generating retroelement (DGR)",
#      the same words the text uses. The figure and the manuscript said two
#      different things for one thing.
#
#   3. "Structural / viral assembly" is now "Structural and other phage genes".
#      The category is assigned to every annotated gene that is not defence or
#      AMG, not only to structural ones, and the old name claimed more.
#
#   4. Categories that no figure draws (transcriptional regulator, SOS response,
#      anti-restriction) were removed. relabel_map() now fails loudly if a
#      script asks for one, instead of carrying dead entries.
# =============================================================================

MAP_CATEGORIES <- c(
  # --- context -------------------------------------------------------------
  "Host (MAG / flanking)",
  "Structural and other phage genes",
  "Hypothetical / unannotated",

  # --- lifecycle accents ---------------------------------------------------
  "Integration machinery",
  "Diversity-generating retroelement (DGR)",

  # --- defence, categories as in Table S4 and Data S2 ----------------------
  "Defence: Nucleic acid degradation",
  "Defence: Genome modification",
  "Defence: Toxin-antitoxin (Abortive infection)",
  "Defence: Other / uncharacterised",

  # --- auxiliary metabolic genes, KEGG categories as in Data S2 ------------
  "AMG: Metabolism of cofactors and vitamins",
  "AMG: Carbohydrate metabolism",
  "AMG: Folding, sorting and degradation",
  "AMG: Energy metabolism"
)

MAP_PALETTE <- c(
  "Host (MAG / flanking)"                          = "#333333",
  "Structural and other phage genes"               = "#4477AA",
  "Hypothetical / unannotated"                     = "#F2F2F2",

  "Integration machinery"                          = "#EE3377",
  "Diversity-generating retroelement (DGR)"        = "#66CCEE",

  "Defence: Nucleic acid degradation"              = "#882255",
  "Defence: Genome modification"                   = "#EE6677",
  "Defence: Toxin-antitoxin (Abortive infection)"  = "#EE8866",
  "Defence: Other / uncharacterised"               = "#999999",

  "AMG: Metabolism of cofactors and vitamins"      = "#EEDD88",
  "AMG: Carbohydrate metabolism"                   = "#117733",
  "AMG: Folding, sorting and degradation"          = "#CCBB44",
  "AMG: Energy metabolism"                         = "#228833"
)

stopifnot(setequal(MAP_CATEGORIES, names(MAP_PALETTE)))
stopifnot(!any(duplicated(MAP_PALETTE)))

# Old spellings are kept as inputs so the map scripts do not have to change:
# whatever they pass in, it comes out in the current vocabulary.
MAP_RELABEL <- c(
  # --- context -------------------------------------------------------------
  "Host (MAG / Flanking)"                   = "Host (MAG / flanking)",
  "Structural / viral assembly"             = "Structural and other phage genes",
  "Structural / Viral Assembly"             = "Structural and other phage genes",
  "Hypothetical / Unannotated"              = "Hypothetical / unannotated",
  "Integration Machinery"                   = "Integration machinery",
  "Hypermutation (DGR)"                     = "Diversity-generating retroelement (DGR)",

  # --- defence: mechanism-first labels tried in Figure 5 -------------------
  "Restriction-modification (Nucleic acid degradation)" = "Defence: Nucleic acid degradation",
  "DNA methylation (Genome modification)"   = "Defence: Genome modification",
  "Toxin-antitoxin (Abortive infection)"    = "Defence: Toxin-antitoxin (Abortive infection)",

  # --- defence: current names and the American spellings they replace ------
  "Defense: Restriction Modification"       = "Defence: Nucleic acid degradation",
  "Defense: Nucleic acid degradation"       = "Defence: Nucleic acid degradation",
  "Defense: DNA Methylation"                = "Defence: Genome modification",
  "Defense: Genome modification"            = "Defence: Genome modification",
  "Defense: Toxin-Antitoxin"                = "Defence: Toxin-antitoxin (Abortive infection)",
  "Defense: Toxin-antitoxin (Abortive infection)" = "Defence: Toxin-antitoxin (Abortive infection)",
  "Defense: Other / uncharacterised"        = "Defence: Other / uncharacterised",

  # --- AMG shorthands ------------------------------------------------------
  "AVG: Folding, sorting and degradation"   = "AMG: Folding, sorting and degradation",
  "AVG: sulfur relay (Folding, sorting & degradation)" = "AMG: Folding, sorting and degradation",
  "AMG: Cofactors & Vitamins"               = "AMG: Metabolism of cofactors and vitamins",
  "AMG: Carbohydrate & Glycan"              = "AMG: Carbohydrate metabolism",

  # --- bare category names, as they appear in Data S2 ----------------------
  "Nucleic acid degradation"                = "Defence: Nucleic acid degradation",
  "Genome modification"                     = "Defence: Genome modification",
  "Other / uncharacterised systems"         = "Defence: Other / uncharacterised",
  "Metabolism of cofactors and vitamins"    = "AMG: Metabolism of cofactors and vitamins",
  "Carbohydrate metabolism"                 = "AMG: Carbohydrate metabolism",
  "Folding, sorting and degradation"        = "AMG: Folding, sorting and degradation",
  "Energy metabolism"                       = "AMG: Energy metabolism"
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
