# ==============================================================================
# 08_Fig5.R — Figure 5. Co-resident prophages in a Vibrio diabolicus MAG.
# Successor of Fig5v3.R. Ready to run.
#
# CHANGES 23/9
#   * Dam is coloured as nucleic acid degradation, the category Data S2 gives
#     it and the category of the complete RM Type II locus it closes with the
#     adjacent mcrA (L01512).
#   * "mec" is labelled "C40 peptidase", as in Table S6 (KEGG K21140).
#   * The HTH regulator is labelled "HTH reg.", so every row of Table S6 has
#     its label in the figure.
#   * MAPS_DIR: the absolute path lives in one variable now.
#
# WHAT CHANGED FROM Fig5v3.R — TWO THINGS, NEITHER TOUCHING THE DATA
#
#   1. TITLES. They read "vOTU NODE_12" and "vOTU NODE_49", which identify
#      nothing: there is a NODE_12 in most assemblies. Figure 4 gives full
#      contig names, so these now do too, and the two figures carry the same
#      information in the same order:
#
#        <letter> | <lifestyle> Caudoviricetes · <identifier> (<size> · <zone>
#                   · <fraction> · <origin>)
#
#      The fraction is added because it is the variable Figure 3e is built on:
#      these two prophages are CUL, the top of that gradient at 56.2% of loci
#      on the exclusion axis, against the BAC and VIR genomes of Figure 4 at
#      29.2% and 10.3%. The origin (culture-dependent) is kept as well: origin
#      and fraction are different variables and neither replaces the other.
#
#      The size shown is the PROPHAGE SPAN, the same number as the bracket
#      drawn inside the panel. The earlier "46.7 kb representative genome" was
#      removed for putting two different kb figures in one panel; that reason
#      still holds, so the vOTU length stays in Table S3 and in the text.
#
#   2. get_legend(). It matched a grob called "guide-box", which ggplot2 3.5
#      renamed to "guide-box-bottom". On a current install `which()` returns
#      integer(0) and the script dies on the index rather than telling you why.
#
# EVERYTHING ELSE IS UNTOUCHED, and deliberately so: the boundaries, the
# offsets, the manually injected ORF and every per-gene category are curated
# against three independent sources and documented gene by gene below. They
# must not be re-derived from a pipeline table, because VIBRANT's coordinates
# are local to the extracted prophage and Bakta's are on the host contig.
#
# ------------------------------------------------------------------------------
# TRAZABILIDAD DE LA EVIDENCIA (todo verificado gen a gen, tres fuentes):
#
#   Bakta gff3   -> esqueleto de coordenadas en el MAG completo (contexto flanqueante)
#   VIBRANT tsv  -> identidad funcional especifica (KO / Pfam / VOG + e-value)
#                   OJO: coordenadas LOCALES del profago extraido. Offsets resueltos
#                   con doble ancla (inicio Y fin coincidentes, error 0 pb):
#                     NODE_12: offset = +62796  (ancla mec, K21140, e=3.5e-16)
#                     NODE_49: offset = +35088  (ancla yoeB, K19158, e=2.2e-45)
#   PHASTEST txt -> sitios att, integrasas y limite de region, de forma independiente
#
# LIMITES: se usa para cada profago la evidencia que realmente es solida.
#
#   NODE_12 -> attL/attR de PHASTEST (58,524-58,558 y 116,432-116,466).
#              Repeticion directa de 35 pb, y attL solapa el extremo 3' del
#              tRNA-Ser(gga) (_00055, 58,467-58,554) = firma clasica de integracion
#              sitio-especifica en un tRNA. El span (57.9 kb) coincide exactamente
#              con el REGION_LENGTH del summary de PHASTEST.
#
#   NODE_49 -> NO se usan los att de PHASTEST. Su attL (13,411-13,422) es un match
#              de solo 12 pb, baja complejidad (CAAAGCGAAAAA), y cae DENTRO de un
#              operon bacteriano. Se usa el limite del vOTU curado
#              (VIBRANT/CheckV) = 35,089-74,498 = 39,410 pb, identico a Table S3.
#
# QUEDA FUERA A PROPOSITO (y ambas herramientas coinciden en excluirlo):
#   NODE_49: bloque RelE/ParD (_01348/_01349) + 2a integrasa (_01350). La hebra
#            se invierte y el bloque lleva integrasa propia -> elemento movil
#            adyacente en tandem. Se dibuja en gris pero SI se etiqueta.
#
# CORRECCIONES DE VERSIONES ANTERIORES:
#   (1) _00069 ya NO se etiqueta como mcrA: solo era la anotacion generica de Bakta.
#   (2) Las dos copias REALES de mcrA (K07451) estan en 73,830-74,198 y 82,710-83,069.
#   (3) La segunda mcrA NO existe en el gff3 (Pyrodigal no la predijo); se inyecta.
#   (4) tusA corregido a _01298 (VIBRANT K04085, e=2.5e-10).
#   (5) Las alturas de etiqueta se asignan con un algoritmo anti-solape generico.
# ==============================================================================

library(ggplot2)
library(dplyr)
library(gggenes)
library(cowplot)

# One line to change for the repository: MAPS_DIR is the only absolute path.
MAPS_DIR  <- "C:/Users/gasch/Documents/R/vOTUs-R/MAPAS"
path_gff3 <- file.path(MAPS_DIR, "lisogenicos/MCHA-OR-22-GLI_bin.6_NODE_49_length_92740_cov_10.407175_pro1/MCHA-OR-22-GLI_bin.6/MCHA-OR-22-GLI_bin.6.gff3")

NAME_12 <- "MCHA-OR-22-GLI_bin.6_NODE_12"
NAME_49 <- "MCHA-OR-22-GLI_bin.6_NODE_49"

# ------------------------------------------------------------------------------
# 1. LIMITES
# ------------------------------------------------------------------------------
ATT_12_L <- 58524
ATT_12_R <- 116466
P12_START <- ATT_12_L
P12_END   <- ATT_12_R

P49_START <- 35089
P49_END   <- 74498

SPAN_12_KB <- round((P12_END - P12_START) / 1000, 1)   # 57.9 kb
SPAN_49_KB <- round((P49_END - P49_START) / 1000, 1)   # 39.4 kb

VIEW_12 <- c(56000, 120000)
VIEW_49 <- c(30000,  80000)

# ------------------------------------------------------------------------------
# 2. CATEGORIAS Y PALETA -- desde palette_maps.R, compartidas por Fig 4, 5 y 6.
# ------------------------------------------------------------------------------
# NO definir aqui categorias_globales ni paleta_global. Antes cada figura tenia
# las suyas y dos colores acabaron significando cosas distintas entre figuras:
#   #ff7f00  AMG cofactores (Fig 4)      vs  toxina-antitoxina (Fig 5 y 6)
#   #4daf4a  AMG carbohidratos (Fig 4)   vs  "Innate Immune System" (Fig 6)
# Las case_when() de mas abajo siguen escribiendo los nombres antiguos:
# relabel_map() los traduce.
source("palette_maps.R")
categorias_globales <- MAP_CATEGORIES
paleta_global       <- MAP_PALETTE

tema_isme <- theme_minimal() + theme(
  panel.grid = element_blank(), axis.line.x = element_line(color = "black"),
  axis.text.y = element_blank(), axis.title.y = element_blank(),
  axis.title.x = element_text(size = 11),
  legend.position = "none", plot.title = element_text(face = "bold", size = 15),
  plot.margin = margin(t = 10, r = 20, b = 15, l = 20)
)

# Same shape as ttl() in 07_Fig4.R, so a reader moving between the two figures
# meets the same fields in the same order.
ttl5 <- function(letter, votu, span_kb) {
  bquote(bold(.(paste0(letter, " | Temperate "))) *
           bolditalic("Caudoviricetes") *
           bold(.(paste0(" prophage \u00b7 ", votu, " in "))) *
           bolditalic("Vibrio diabolicus") *
           bold(.(paste0(" MAG (", span_kb, " kb \u00b7 OR \u00b7 CUL \u00b7 culture-dependent)"))))
}

# ------------------------------------------------------------------------------
# 3. PARSER + INYECCION DE ORFs QUE BAKTA NO PREDIJO
# ------------------------------------------------------------------------------
parse_gff3 <- function(file, target_contig) {
  lines <- readLines(file)
  fasta_idx <- grep("^##FASTA", lines)
  if (length(fasta_idx) > 0) lines <- lines[1:(fasta_idx[1] - 1)]
  lines <- lines[!grepl("^#", lines) & trimws(lines) != ""]
  df <- read.table(text = lines, sep = "\t", quote = "", stringsAsFactors = FALSE,
                   fill = TRUE, flush = TRUE)[, 1:9]
  df %>%
    filter(V1 == target_contig, V3 %in% c("CDS", "tRNA")) %>%
    mutate(
      locus = sub(".*locus_tag=([^;]+).*", "\\1", V9),
      product = ifelse(grepl("product=", V9),
                       sub(".*product=([^;]+).*", "\\1", V9), "hypothetical protein"),
      start = as.numeric(V4), end = as.numeric(V5),
      strand = ifelse(V7 == "-", -1, 1)
    ) %>%
    select(locus, product, start, end, strand)
}

# ORFs confirmados por VIBRANT + PHASTEST que Pyrodigal/Bakta NO llamo.
orfs_rescatados_c1 <- data.frame(
  locus   = "MANUAL_mcrA2",
  product = "mcrA (5-methylcytosine-specific restriction enzyme A); VIBRANT K07451 e=2.3e-17; PHASTEST e=5.18e-37; ORF no predicho por Bakta",
  start   = 82710, end = 83069, strand = 1,
  stringsAsFactors = FALSE
)

preparar <- function(df) {
  df %>% arrange(start) %>%
    mutate(direccion_hebra = strand == 1, midpoint = (start + end) / 2)
}

# ------------------------------------------------------------------------------
# 4. ALGORITMO ANTI-SOLAPE DE ETIQUETAS
# ------------------------------------------------------------------------------
asignar_niveles <- function(df, min_dist, niveles = c(1.17, 1.27, 1.37, 1.47)) {
  df$label_y <- NA_real_
  idx <- which(df$label != "")
  if (!length(idx)) return(df)
  idx <- idx[order(df$start[idx])]
  ultimo <- rep(-Inf, length(niveles))
  for (i in idx) {
    libres <- which(df$start[i] - ultimo >= min_dist)
    k <- if (length(libres)) min(libres) else which.min(ultimo)
    df$label_y[i] <- niveles[k]
    ultimo[k] <- df$start[i]
  }
  df
}

# ------------------------------------------------------------------------------
# 5. PANEL A -- NODE_12 (contig_1)
# ------------------------------------------------------------------------------
mapa_12 <- bind_rows(parse_gff3(path_gff3, "contig_1"), orfs_rescatados_c1) %>%
  preparar() %>%
  mutate(
    tag = sub(".*_(\\d{5})$", "\\1", locus),
    Category = relabel_map(case_when(
      locus == "MANUAL_mcrA2"             ~ "Defense: Restriction Modification",
      tag %in% c("00056","00131","00133") ~ "Integration Machinery",
      tag == "00064"                      ~ "Defense: Phage Communication",   # luxR-like, ver nota
      # Dam takes the category Data S2 gives it, "Nucleic acid degradation"
      # (PADLOC MTase_II), which is also the category of the locus it closes:
      # with the mcrA next to it, L01512 is a COMPLETE, confirmed RM Type II.
      tag == "00080"                      ~ "Defense: Restriction Modification", # Dam, PF05869.11, e=2.3e-37
      tag == "00082"                      ~ "Defense: Restriction Modification", # mcrA, K07451, e=5e-11
      tag == "00116"                      ~ "AMG: Folding, sorting and degradation", # mec, K21140, e=3.5e-16
      tag %in% c("00128","00129")         ~ "Hypermutation (DGR)",             # Avd / RT
      grepl("hypothetical", product, ignore.case = TRUE) ~ "Hypothetical / Unannotated",
      midpoint < P12_START | midpoint > P12_END          ~ "Host (MAG / Flanking)",
      TRUE                                               ~ "Structural / Viral Assembly"
    )),
    
    label = case_when(
      tag == "00055"          ~ "tRNA-Ser",   # attB
      tag == "00056"          ~ "int",
      # _00057 y _00074 (transposasas) YA NO se etiquetan: estaban rotuladas y
      # definidas en la leyenda pero no se discuten en ningun punto del texto.
      # Un gen etiquetado y no discutido es una pregunta gratis para el revisor.
      # _00064 YA NO se etiqueta, por la misma razon. Estaba rotulado "luxR" en
      # cursiva, lo que afirma un nombre de gen concreto. No lo tiene: Bakta lo
      # llama "helix-turn-helix transcriptional regulator" sin gene name, el KO
      # K03556 (malT, 901 aa) es imposible en un ORF de 141 aa, y blastp contra
      # nr devuelve "response regulator transcription factor" de Vibrio al 97%
      # de identidad sobre 141 aa. Ni malT, ni LuxR de quorum sensing, ni
      # arbitrium. El color sigue mostrando que ahi hay un regulador; la
      # etiqueta afirmaba mas de lo que se puede sostener.
      # Labelled as the neutral family name documented in Table S6, not as a
      # gene name: Bakta gives it no gene name and blastp returns a Vibrio
      # response regulator. The old "luxR" label asserted a system.
      tag == "00064"          ~ "HTH reg.",
      tag == "00080"          ~ "Dam",
      tag == "00082"          ~ "mcrA",
      locus == "MANUAL_mcrA2" ~ "mcrA",
      tag == "00116"          ~ "C40 peptidase",   # K21140; Table S6
      tag == "00128"          ~ "Avd",
      tag == "00129"          ~ "RT",
      tag == "00131"          ~ "int",
      tag == "00133"          ~ "rec",
      TRUE ~ ""
    )
  ) %>%
  asignar_niveles(min_dist = 3600)

pA <- ggplot(mapa_12, aes(xmin = start, xmax = end, y = 1,
                          fill = Category, forward = direccion_hebra)) +
  annotate("rect", xmin = P12_START, xmax = P12_END, ymin = 0.86, ymax = 1.13,
           fill = "grey93", color = NA) +
  geom_gene_arrow(arrow_body_height = unit(5, "mm"), arrowhead_height = unit(6, "mm"),
                  linewidth = 0.2, color = "black") +
  geom_segment(data = subset(mapa_12, label != ""),
               aes(x = midpoint, xend = midpoint, y = 1.09, yend = label_y - 0.012),
               inherit.aes = FALSE, color = "grey55", linewidth = 0.25) +
  geom_text(data = subset(mapa_12, label != ""),
            aes(x = midpoint, y = label_y, label = label),
            inherit.aes = FALSE, hjust = 0.5, vjust = 0, size = 3.9, fontface = "italic") +
  annotate("segment", x = P12_START, xend = P12_START, y = 0.72, yend = 1.05,
           color = "grey40", linewidth = 0.8, linetype = "dotted") +
  annotate("segment", x = P12_END, xend = P12_END, y = 0.72, yend = 1.05,
           color = "grey40", linewidth = 0.8, linetype = "dotted") +
  annotate("point", x = c(P12_START, P12_END), y = 0.80, shape = 25, size = 2.6,
           fill = "black", color = "black") +
  annotate("text", x = c(P12_START + 900, P12_END + 900), y = 0.755,
           label = c("attL", "attR"), size = 2.9, fontface = "italic") +
  annotate("segment", x = P12_START, xend = P12_END, y = 1.62, yend = 1.62, linewidth = 0.5) +
  annotate("segment", x = c(P12_START, P12_END), xend = c(P12_START, P12_END),
           y = 1.58, yend = 1.66, linewidth = 0.5) +
  # Etiqueta del corchete NEUTRA e IDENTICA en ambos paneles. El criterio usado
  # en cada uno (att vs. limite curado) se explica en la leyenda y en M&M.
  annotate("text", x = mean(c(P12_START, P12_END)), y = 1.70,
           label = paste0("Prophage region: ", SPAN_12_KB, " kb"),
           size = 3.1, fontface = "italic") +
  scale_x_continuous(breaks = seq(55000, 120000, 5000),
                     labels = function(x) paste0(x / 1000, " kb")) +
  coord_cartesian(xlim = VIEW_12, ylim = c(0.66, 1.80)) +
  scale_fill_manual(values = paleta_global, drop = FALSE) +
  tema_isme +
  labs(x = "Position in host MAG (kb)",
       title = ttl5("a", NAME_12, SPAN_12_KB))

# ------------------------------------------------------------------------------
# 6. PANEL B -- NODE_49 (contig_13)
# ------------------------------------------------------------------------------
mapa_49 <- parse_gff3(path_gff3, "contig_13") %>%
  preparar() %>%
  mutate(
    tag = sub(".*_(\\d{5})$", "\\1", locus),
    Category = relabel_map(case_when(
      tag == "01323"                      ~ "Integration Machinery",   # int, K14059, e=1.6e-18
      tag %in% c("01292","01321","01322") ~ "Defense: Toxin-Antitoxin", # chpB / yoeB / yefM
      tag == "01298"                      ~ "AMG: Folding, sorting and degradation", # tusA, K04085
      tag %in% c("01303","01304")         ~ "Hypermutation (DGR)",      # RT / Avd
      grepl("hypothetical", product, ignore.case = TRUE) ~ "Hypothetical / Unannotated",
      midpoint < P49_START | midpoint > P49_END          ~ "Host (MAG / Flanking)",
      TRUE                                               ~ "Structural / Viral Assembly"
    )),
    
    label = case_when(
      # _01290 (20,805-34,643, ORF de 13.8 kb) NO se etiqueta a proposito: su
      # identidad esta en disputa entre las tres herramientas Y esta fuera del
      # limite. _01356 (norM) tampoco: gen de eflujo del hospedador, 7 kb fuera.
      tag == "01292" ~ "chpB",
      tag == "01298" ~ "tusA",
      tag == "01301" ~ "tail",
      tag == "01303" ~ "RT",
      tag == "01304" ~ "Avd",
      tag == "01321" ~ "yoeB",
      tag == "01322" ~ "yefM",
      tag == "01323" ~ "int",
      # Estos tres SI se etiquetan aunque esten fuera: identidad inequivoca y
      # relevancia directa para la historia de lisogenia.
      tag == "01348" ~ "RelE",
      tag == "01349" ~ "ParD",
      tag == "01350" ~ "int*",
      TRUE ~ ""
    )
  ) %>%
  asignar_niveles(min_dist = 3000)

MGE_START <- 74514   # _01346, primer gen tras el limite curado
MGE_END   <- 76955   # _01350, 2a integrasa

pB <- ggplot(mapa_49, aes(xmin = start, xmax = end, y = 1,
                          fill = Category, forward = direccion_hebra)) +
  annotate("rect", xmin = P49_START, xmax = P49_END, ymin = 0.86, ymax = 1.13,
           fill = "grey93", color = NA) +
  geom_gene_arrow(arrow_body_height = unit(5, "mm"), arrowhead_height = unit(6, "mm"),
                  linewidth = 0.2, color = "black") +
  geom_segment(data = subset(mapa_49, label != ""),
               aes(x = midpoint, xend = midpoint, y = 1.09, yend = label_y - 0.012),
               inherit.aes = FALSE, color = "grey55", linewidth = 0.25) +
  geom_text(data = subset(mapa_49, label != ""),
            aes(x = midpoint, y = label_y, label = label),
            inherit.aes = FALSE, hjust = 0.5, vjust = 0, size = 3.9, fontface = "italic") +
  annotate("segment", x = P49_START, xend = P49_START, y = 0.72, yend = 1.05,
           color = "grey40", linewidth = 0.8, linetype = "dotted") +
  annotate("segment", x = P49_END, xend = P49_END, y = 0.72, yend = 1.05,
           color = "grey40", linewidth = 0.8, linetype = "dotted") +
  annotate("segment", x = P49_START, xend = P49_END, y = 1.62, yend = 1.62, linewidth = 0.5) +
  annotate("segment", x = c(P49_START, P49_END), xend = c(P49_START, P49_END),
           y = 1.58, yend = 1.66, linewidth = 0.5) +
  annotate("text", x = mean(c(P49_START, P49_END)), y = 1.70,
           label = paste0("Prophage region: ", SPAN_49_KB, " kb"),
           size = 3.1, fontface = "italic") +
  # Corchete secundario, mas fino, para el elemento adyacente: lo situa sin
  # sugerir que forme parte del profago.
  annotate("segment", x = MGE_START, xend = MGE_END, y = 0.82, yend = 0.82,
           linewidth = 0.4, color = "grey30") +
  annotate("segment", x = c(MGE_START, MGE_END), xend = c(MGE_START, MGE_END),
           y = 0.79, yend = 0.85, linewidth = 0.4, color = "grey30") +
  annotate("text", x = mean(c(MGE_START, MGE_END)) + 900, y = 0.735,
           label = "adjacent MGE", size = 2.9, fontface = "italic", color = "grey30") +
  scale_x_continuous(breaks = seq(30000, 80000, 5000),
                     labels = function(x) paste0(x / 1000, " kb")) +
  coord_cartesian(xlim = VIEW_49, ylim = c(0.66, 1.80)) +
  scale_fill_manual(values = paleta_global, drop = FALSE) +
  tema_isme +
  labs(x = "Position in host MAG (kb)",
       title = ttl5("b", NAME_49, SPAN_49_KB))

# ------------------------------------------------------------------------------
# 7. LEYENDA
# ------------------------------------------------------------------------------
p_legend <- ggplot(bind_rows(mapa_12, mapa_49),
                   aes(xmin = start, xmax = end, y = 1, fill = Category)) +
  geom_gene_arrow() + scale_fill_manual(values = paleta_global, drop = FALSE) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold"),
        legend.key.size = unit(0.7, "cm"), legend.spacing.x = unit(0.5, "cm"),
        legend.box.margin = margin(t = 15, b = 15)) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE,
                             keywidth = unit(0.9, "cm"), keyheight = unit(0.9, "cm")))

# ggplot2 3.5 renamed the legend grob from "guide-box" to "guide-box-bottom".
# The previous version indexed g$grobs by that name, so on a current install
# which() returns integer(0) and the script dies on the subscript instead of
# saying what is wrong. Current accessor first, old one as fallback.
get_legend <- function(p) {
  g <- try(cowplot::get_plot_component(p, "guide-box-bottom", return_all = TRUE),
           silent = TRUE)
  if (!inherits(g, "try-error") && !is.null(g) && !inherits(g, "zeroGrob")) return(g)
  cowplot::get_legend(p)
}
legend_final <- get_legend(p_legend)
stopifnot(!is.null(legend_final))

figura_completa <- plot_grid(pA, pB, legend_final, ncol = 1, rel_heights = c(1, 1, 0.55))

# ------------------------------------------------------------------------------
# 8. EXPORTACION + CHEQUEO EN CONSOLA
# ------------------------------------------------------------------------------
cat("\n--- Etiquetas panel a (NODE_12) ---\n")
print(mapa_12 %>% filter(label != "") %>%
        transmute(label, start, end, strand,
                  dentro = start >= P12_START & end <= P12_END, Category))
cat("\n--- Etiquetas panel b (NODE_49) ---\n")
print(mapa_49 %>% filter(label != "") %>%
        transmute(label, start, end, strand,
                  dentro = start >= P49_START & end <= P49_END, Category))
cat("\n--- Categorias dibujadas (deben coincidir con la leyenda de la Fig. 4) ---\n")
print(bind_rows(mapa_12, mapa_49) %>% count(Category, .drop = FALSE) %>% filter(n > 0))

# ------------------------------------------------------------------------------
# EXPORT — PNG for looking at, TIFF and EPS for the journal
#
# EPS via cairo_ps rather than the base postscript device: the base one cannot
# embed fonts or handle transparency and silently drops what it cannot draw.
# cairo_ps is a build-time option in R, so it is checked rather than assumed;
# on a build without it the EPS is skipped with a message instead of erroring
# after the other two files are already on disk.
#
# This figure uses no alpha (the prophage band is a solid grey93), so nothing
# should rasterise. fallback_resolution is set anyway so that adding a
# semi-transparent layer later degrades to a 600-dpi bitmap in that region
# rather than disappearing.
# ------------------------------------------------------------------------------
save_all <- function(plot, base, w, h) {
  ggsave(paste0(base, ".png"),  plot, width = w, height = h, dpi = 600, bg = "white")
  ggsave(paste0(base, ".tiff"), plot, width = w, height = h, dpi = 600,
         compression = "lzw", bg = "white")
  if (isTRUE(capabilities("cairo"))) {
    ggsave(paste0(base, ".eps"), plot, width = w, height = h,
           device = cairo_ps, fallback_resolution = 600, bg = "white")
    cat("  wrote ", base, ".png / .tiff / .eps\n", sep = "")
  } else {
    warning("This R build has no cairo support, so no EPS was written. ",
            "Install R with cairo, or open the TIFF in Illustrator/Inkscape ",
            "and export EPS from there.", call. = FALSE)
    cat("  wrote ", base, ".png / .tiff  (EPS skipped: no cairo)\n", sep = "")
  }
}

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
save_all(figura_completa, paste0("Figure_5_LysogenicMaps_", stamp), 16, 9.5)
save_all(figura_completa, "Figure_5_LysogenicMaps", 16, 9.5)
cat("\nSend THIS file: Figure_5_LysogenicMaps_", stamp, ".png\n", sep = "")
