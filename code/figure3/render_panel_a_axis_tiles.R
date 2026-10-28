#!/usr/bin/env Rscript
# Render Figure 3 panel A tiles: the full hub->TLS axis (Hub-only / Hub+TLS
# transition / E-TLS / PFL-TLS / SFL-TLS) for two public samples (labels from the
# main pipeline outputs) and two in-house Visium sections (classified here with
# the exact rules of spatial/run_full_tls_pipeline.R).

library(Seurat); library(ggplot2); library(grid); library(dplyr)

OUT_DIR <- "data/hub_tls_paper_outline/figure3_assembled/panel_a_tiles"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

STAGE_COLORS <- c("E-TLS" = "#D24B40", "PFL-TLS" = "#E7A03D", "SFL-TLS" = "#5DA05D")
# discrete spot labelling (user preference over the continuous hub field):
# Hub-only (blue) / Hub+TLS overlap (purple) / TLS (red). The overlap is kept
# visible because ~40% of classified spots are double-positive and it reads as
# the transition zone between the two poles.
MAP_COLORS <- c("Hub" = "#6C8EBF", "Hub+TLS" = "#7A5195", "TLS" = "#C41E3A")

render_tile <- function(seu, tls_stage, is_hub, out_png) {
  coords <- GetTissueCoordinates(seu)
  cls0 <- ifelse(is_hub & !is.na(tls_stage), "Hub+TLS",
          ifelse(is_hub, "Hub",
          ifelse(!is.na(tls_stage), "TLS", NA)))
  cls0 <- setNames(cls0, colnames(seu))
  keep <- intersect(rownames(coords), colnames(seu))   # coords may drop off-tissue spots
  # all spots as light-grey lattice (tissue outline), classified spots on top
  # with a white halo so they lift off the lattice
  df <- data.frame(x = coords[keep, 1], y = coords[keep, 2], cls = unname(cls0[keep]),
                   stringsAsFactors = FALSE)
  bg <- df[, c("x", "y")]
  fg <- df[!is.na(df$cls), ]
  fg$cls <- factor(fg$cls, levels = names(MAP_COLORS))
  fg <- fg[order(fg$cls), ]              # Hub (bottom) -> Hub+TLS -> TLS (top)
  # point sizes matched to the Visium lattice: diameter ~= 1.1x NN spacing
  sub <- bg[sample(nrow(bg), min(1500, nrow(bg))), ]
  d <- as.matrix(dist(sub[, c("x", "y")])); diag(d) <- NA
  spacing <- median(apply(d, 2, min, na.rm = TRUE))
  sz_bg <- 0.55 * spacing / diff(range(bg$x)) * 98   # lattice dots nearly touching
  sz_fg <- 1.05 * spacing / diff(range(bg$x)) * 98
  bar <- 5 * spacing
  bx0 <- min(bg$x) + 0.02 * diff(range(bg$x)); by <- max(bg$y) + 0.06 * diff(range(bg$y))
  p <- ggplot() +
    geom_point(data = bg, aes(x, y), colour = "grey90", size = sz_bg) +
    geom_point(data = fg, aes(x, y), colour = "white", size = sz_fg * 1.35, stroke = 0) +
    geom_point(data = fg, aes(x, y, colour = cls), size = sz_fg) +
    annotate("segment", x = bx0, xend = bx0 + bar, y = by, yend = by,
             colour = "black", linewidth = 1.2) +
    annotate("text", x = bx0 + bar / 2, y = by, label = "500 µm", vjust = -0.7,
             size = 7 / .pt, family = "Arial") +
    scale_colour_manual(values = MAP_COLORS, guide = "none") +
    coord_fixed() + theme_void() +
    theme(plot.margin = margin(2, 2, 2, 2, "mm"))
  ggsave(out_png, p, width = 98, height = 98, units = "mm", dpi = 300, bg = "white")
  table(fg$cls)
}

# public samples: labels already computed by the main pipeline
## Chosen by screen_panel_a_samples.R across all 158 samples: highest TLS
## clustering + hub-adjacency-to-TLS (adj_score 2.29 / 1.73). Both are
## BreastCancer_Spatial_CID (Visium); the previous PDAC/Lung pair were confetti.
for (spec in list(
    list(rds = "data/spatial/tls_results_v3/BreastCancer_Spatial_CID/CID4535_tls.Rds",
         tag = "Breast_CID4535"),
    list(rds = "data/spatial/tls_results_v3/BreastCancer_Spatial_CID/CID4465_tls.Rds",
         tag = "Breast_CID4465"))) {
  seu <- readRDS(spec$rds)
  tab <- render_tile(seu, seu$tls_stage, seu$is_immune_hub,
                     file.path(OUT_DIR, paste0(spec$tag, "__hub_tls_axis.png")))
  cat(spec$tag, ":", paste(names(tab), tab, sep = "=", collapse = " "), "\n")
}

# in-house samples: classify with the pipeline rules
detection_modules <- list(
  Meylan_5    = c("CXCL13", "CCL19", "CCL21", "LAMP3", "CCR7"),
  Cabrita_12  = c("CD3D", "CD3E", "CD8A", "CD19", "CD79A", "MS4A1",
                  "CXCL13", "CCL19", "CCL21", "LAMP3", "CCR7", "ICOS"),
  Sautes_9    = c("CD3D", "CD3E", "CD8A", "CD19", "CD79A", "MS4A1",
                  "CXCL13", "CCL19", "CCL21"),
  Bcell_mod   = c("CD19", "CD79A", "MS4A1", "BCL6", "AICDA", "CD10",
                  "MKI67", "CR2", "CD21", "CD23", "FDCSP"),
  Tcell_mod   = c("CD3D", "CD3E", "CD8A", "ICOS", "CXCR5", "PDCD1",
                  "CD40LG", "PRDM1", "IRF4"),
  Chemo_core  = c("CXCL13", "CCL19", "CCL21", "MADCAM1", "ICAM1", "VCAM1"))
maturation_modules <- list(
  E_TLS_GENES   = c("CCL19", "CCL21", "CXCL13", "LAMP3", "CCR7"),
  PFL_TLS_GENES = c("CD3D", "CD3E", "CD8A", "CD19", "CD79A", "MS4A1",
                    "ICOS", "CXCR5", "PDCD1", "CD40LG"),
  SFL_TLS_GENES = c("BCL6", "AICDA", "CD10", "MKI67", "CR2", "CD21",
                    "CD23", "FDCSP", "MZB1", "JCHAIN"))
hub_modules <- list(Hub_chemo = c("CXCL10", "CXCL11", "CXCL9"),
                    Hub_tcell = c("IFNG", "CD3E", "CD8A", "CXCL13"))
stage_map <- c(E_TLS_GENES = "E-TLS", PFL_TLS_GENES = "PFL-TLS", SFL_TLS_GENES = "SFL-TLS")

score_modules <- function(seu, mod_list) {
  for (mname in names(mod_list)) {
    g <- intersect(mod_list[[mname]], rownames(seu))
    seu@meta.data[[mname]] <- if (length(g) == 0) 0 else
      AddModuleScore(seu, features = list(g), name = mname,
                     ctrl = min(50, max(5, length(g) * 5)))@meta.data[[paste0(mname, "1")]]
  }
  seu
}

for (spec in list(list(rds = "data/ai/Re.rds", tag = "Own_Re"),
                  list(rds = "data/ai/NoRe.rds", tag = "Own_NoRe"))) {
  seu <- readRDS(spec$rds)
  if (!"data" %in% slotNames(seu@assays[[DefaultAssay(seu)]]) ||
      identical(GetAssayData(seu, layer = "data"), GetAssayData(seu, layer = "counts")))
    seu <- NormalizeData(seu, verbose = FALSE)
  seu <- score_modules(seu, detection_modules)
  seu <- score_modules(seu, maturation_modules)
  seu <- score_modules(seu, hub_modules)
  seu$composite_score <- rowMeans(as.matrix(seu@meta.data[, names(detection_modules)]), na.rm = TRUE)
  seu$tls_spot <- seu$composite_score >= quantile(seu$composite_score, 0.90, na.rm = TRUE)
  seu$hub_score <- rowMeans(as.matrix(seu@meta.data[, c("Hub_chemo", "Hub_tcell")]), na.rm = TRUE)
  seu$is_hub <- seu$hub_score >= quantile(seu$hub_score, 0.85, na.rm = TRUE)
  smat <- as.matrix(seu@meta.data[, c("E_TLS_GENES", "PFL_TLS_GENES", "SFL_TLS_GENES")])
  assign_v <- stage_map[colnames(smat)[max.col(smat, ties.method = "first")]]
  seu$tls_stage <- NA_character_
  seu$tls_stage[seu$tls_spot] <- assign_v[seu$tls_spot]
  tab <- render_tile(seu, seu$tls_stage, seu$is_hub,
                     file.path(OUT_DIR, paste0(spec$tag, "__hub_tls_axis.png")))
  cat(spec$tag, ":", paste(names(tab), tab, sep = "=", collapse = " "), "\n")
}
