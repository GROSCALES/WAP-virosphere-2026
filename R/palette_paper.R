# =============================================================================
# palette_paper.R — ONE palette for the whole paper. Sourced by every figure
# script; no figure defines its own colours.
#
# THE RULE: ONE HUE FAMILY PER VARIABLE. A reader must be able to learn the
# colour code once and carry it across every panel of every figure.
#
#   green / gold   WHERE the sequence came from   (Origin, and Fraction
#                                                  nested inside it)
#   purple         WHAT the phage does            (Lifestyle)
#   red            WHICH DEFENCE AXIS             (Self-protection vs Exclusion)
#   grey           structural grouping            (Class: Defence vs AMGs)
#   blue           continuous abundance           (heatmaps)
#
# WHY THIS FILE EXISTS — two real problems it fixes:
#
#   1. THE SAME GREEN MEANT THREE DIFFERENT THINGS. Within Figure 2, green
#      was GLI/LIS (fraction) in panels a-b, Culture-dependent (origin) in
#      panel c, and Temperate (lifestyle) in panel d. Three variables sharing
#      one hue, with all three legends assembled onto the same page. Lifestyle
#      is therefore moved to purple, which no other variable uses.
#
#   2. THE TWO FIGURES DISAGREED WITH EACH OTHER. Culture-dependent was
#      #229954 in Figure 2 and #3B7A57 in Figure 3; Temperate was #117864 in
#      Figure 2 and #82B990 in Figure 3 -- the same concept in different
#      colours depending on which figure you were looking at. Defining the
#      palette in one place makes that class of drift impossible.
#
# Fraction deliberately stays inside the Origin hues: CB22/VIRFC are golds
# because they ARE culture-independent, GLI/LIS are greens because they ARE
# culture-dependent. That nesting is informative and is kept.
# =============================================================================

# --- Origin: where the sequence came from -----------------------------------
col_ori <- c(
  "Culture-dependent"   = "#1B7A4B",
  "Culture-independent" = "#C9A227"
)

# --- Fraction: nested inside Origin, same hues, lighter/darker ---------------
col_frac <- c(
  "CB22"  = "#B7950B",   # gold, culture-independent
  "VIRFC" = "#F1C40F",   # gold, culture-independent
  "GLI"   = "#196F3D",   # green, culture-dependent
  "LIS"   = "#58D68D"    # green, culture-dependent
)

# --- Lifestyle: its own hue family, shared with nothing ----------------------
col_life <- c(
  "Temperate" = "#54278F",
  "Virulent"  = "#BCBDDC"
)

# --- Defence axis -----------------------------------------------------------
# NAMING. The axes used to be called "Counter-defence" and "Defence", which
# collided with "Defence" as the name of the gene CLASS (Defence vs AMG). The
# same word then meant two different things at two different levels: a locus
# on the "Defence axis" was a subset of the "Defence class", and a locus on
# the "Counter-defence axis" was ALSO in the Defence class. Renamed to what
# each axis actually does, so the class keeps its name and the axes stop
# borrowing it:
#     Self-protection  the phage shields its own genome from host systems
#                      (orphan methyltransferases; mapa_v5 axis 1)
#     Exclusion        the phage acts against other incoming genetic material
#                      (nucleases, RM, retrons, TA/Abi; mapa_v5 axis 2)
#
# COLOUR. The first draft used purple for axis 1, which is the exact failure
# this file exists to prevent: purple already means Temperate. Red is the only
# family no other variable uses (green/gold = origin, purple = lifestyle,
# blue = heatmaps, grey = class). Light/dark within one hue, the same
# construction as col_life, so the pair reads as two levels of one variable.
col_axis <- c(
  "Self-protection" = "#F3C9C2",   # pale brick — the dominant mass
  "Exclusion"       = "#B03A2E"    # brick     — the rare, informative one
)

# --- Class (Defence vs AMGs): structural, so deliberately neutral -----------
# Previously #5D8596 / #887094 -- the second is a purple-grey and would now
# compete with Lifestyle. Greys carry the grouping without claiming a hue.
# Both spellings are provided: the code has said AVG for a year, the
# manuscript says AMG, and a figure that looks up the wrong key gets grey NA
# bars with no error. Cheaper to carry two names than to debug that once.
col_class <- c(
  "Defence" = "#4D4D4D",
  "AMGs"    = "#BABABA",
  "AMG"     = "#BABABA",
  "AVGs"    = "#BABABA"
)

# --- Fraction display labels, so they cannot drift between figures either ---
frac_lbl <- c("CB22" = "BAC", "GLI" = "CUL", "VIRFC" = "VIR", "LIS" = "LIS")

# --- Fraction display ORDER -------------------------------------------------
# Figure 2 sorts alphabetically (BAC, CUL, LIS, VIR) because there the panels
# are descriptive. Figure 3 sorts by distance from a living cell, because
# there that ordering IS the result: the exclusion axis rises monotonically
# along it. Both orders are named here so neither gets retyped by hand.
frac_lev_alpha <- c("BAC", "CUL", "LIS", "VIR")
frac_lev_cell  <- c("VIR", "LIS", "BAC", "CUL")   # free particles -> cultured cells
