#!/usr/bin/env Rscript
# Render Panel B of Figure 3 (spatial DPT strip) with no titles and enlarged legends.
# Output overwrites panel_a_tiles/spatial_dpt_strip.png so assemble_figure3.R picks
# it up unchanged (the script reuses the cached PNG when present).

library(data.table); library(ggplot2); library(patchwork); library(RANN)
library(ragg); library(png); library(cowplot); library(svglite)

BASE_FONT <- "Arial"
OUT_PNG <- "data/hub_tls_paper_outline/figure3_assembled/panel_a_tiles/spatial_dpt_strip.png"
OUT_SVG <- "data/hub_tls_paper_outline/figures_AI_ready/panel_B_spatial_dpt_standalone.svg"
OUT_PDF <- "data/hub_tls_paper_outline/figures_AI_ready/panel_B_spatial_dpt_standalone.pdf"

DPT_CSV <- "data/spatial/hub_precursor_analysis/paga_dpt_hub_tls_spatial/dpt_pseudotime.csv"
SLING_CSV <- "data/spatial/hub_precursor_analysis/continuous_pseudotime/scheme_a_independent/continuous_pseudotime_results.csv"

CLUSTER_COLORS <- c("Hub_only" = "#43938C",
                    "Hub_TLS_transition" = "#E9B44C",
                    "TLS" = "#C05746")
CLUSTER_LABELS <- c("Hub_only" = "Hub-only",
                    "Hub_TLS_transition" = "Hub-TLS transition",
                    "TLS" = "TLS")

RING_COLORS <- c("0-100"    = "#C05746",
                 "100-300"  = "#E9A13B",
                 "300-500"  = "#5DA05D",
                 "500-1000" = "#43938C",
                 ">1000"    = "#6C8EBF")
RING_LABELS <- names(RING_COLORS)

# Load + derive
dpt_dt <- fread(DPT_CSV)
subcluster_dt <- fread(SLING_CSV, select = c("barcode", "slingshot_cluster"))
dpt_dt <- merge(dpt_dt, subcluster_dt, by = "barcode", all.x = TRUE)
dpt_dt[, cluster := factor(cluster, levels = names(CLUSTER_COLORS))]
dpt_dt[, is_tls_core := slingshot_cluster %in% c("TLS_only", "SFL-TLS")]

dpt_dt[, dist_to_tls_core := {
  core <- .SD[is_tls_core == TRUE, .(coord_x, coord_y)]
  if (nrow(core) == 0) rep(NA_real_, .N)
  else as.vector(RANN::nn2(core, .SD[, .(coord_x, coord_y)], k = 1)$nn.dists)
}, by = sample]

dpt_dt[, dist_ring := cut(dist_to_tls_core,
                          breaks = c(0, 100, 300, 500, 1000, Inf),
                          labels = RING_LABELS, include.lowest = TRUE, right = FALSE)]
dpt_dt[, dist_ring := factor(dist_ring, levels = RING_LABELS)]

# Theme
theme_b <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = BASE_FONT) +
    theme(axis.line = element_line(linewidth = 0.4, colour = "black"),
          axis.ticks = element_line(linewidth = 0.4, colour = "black"),
          axis.text = element_text(size = base_size, colour = "black"),
          axis.title = element_text(size = base_size, colour = "black"),
          plot.title = element_blank(),
          plot.subtitle = element_blank(),
          plot.margin = margin(1, 1, 1, 1, "mm"))
}

# Left: CTC21P spatial cluster map with TLS-core rings
REP_SAMPLE <- "CTC21P"
samp_dt <- dpt_dt[sample == REP_SAMPLE]
core_spots <- samp_dt[is_tls_core == TRUE]
cx <- if (nrow(core_spots) > 0) median(core_spots$coord_x) else median(samp_dt$coord_x)
cy <- if (nrow(core_spots) > 0) median(core_spots$coord_y) else median(samp_dt$coord_y)

circle_coords <- function(cx, cy, r, n = 120) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.table(x = cx + r * cos(theta), y = cy + r * sin(theta), r = as.character(r))
}
ring_circles <- rbindlist(lapply(c(100, 300, 500), function(r) circle_coords(cx, cy, r)))

p_left <- ggplot() +
  geom_point(data = samp_dt,
             aes(x = coord_x, y = coord_y, colour = cluster),
             size = 0.7, alpha = 0.85, stroke = 0) +
  geom_path(data = ring_circles,
            aes(x = x, y = y, group = r),
            colour = "grey20", linewidth = 0.55, linetype = "dashed") +
  scale_colour_manual(values = CLUSTER_COLORS, labels = CLUSTER_LABELS, name = NULL) +
  coord_fixed(expand = FALSE) +
  labs(x = "X (px)", y = "Y (px)") +
  theme_b() +
  theme(axis.text = element_text(size = 6.5),
        axis.title = element_text(size = 6.5),
        legend.position = "none")

# Right: violin of DPT by distance ring
p_right <- ggplot(dpt_dt[!is.na(dist_ring)],
                  aes(x = dist_ring, y = dpt_pseudotime, fill = dist_ring)) +
  geom_violin(scale = "width", alpha = 0.7, linewidth = 0.3, trim = TRUE, colour = NA) +
  geom_boxplot(width = 0.18, outlier.size = 0.2, outlier.alpha = 0.25,
               linewidth = 0.3, fill = "white", colour = "grey20") +
  scale_fill_manual(values = RING_COLORS, name = NULL) +   # title redundant: x-axis already says "Distance ... (px)"
  scale_y_continuous(limits = c(0.15, 0.85), breaks = c(0.2, 0.4, 0.6, 0.8)) +
  labs(x = "Distance to nearest TLS core (px)", y = "Spatial DPT") +
  theme_b() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1, size = 6.5),
        axis.title.y = element_text(size = 6.5, margin = margin(r = 1, unit = "mm")))

# Combine + save
# Let each subplot carry its own legend at the bottom; patchwork collects them.
p_left_legged <- p_left +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.size = unit(3.2, "mm"),
        legend.key.spacing.x = unit(1.5, "mm"),
        legend.text = element_text(size = 6.5, family = BASE_FONT),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0)) +
  guides(colour = guide_legend(override.aes = list(size = 2.8, alpha = 1), nrow = 1))

p_right_legged <- p_right +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 6.5, family = BASE_FONT),
        legend.key.size = unit(3.0, "mm"),
        legend.key.spacing.x = unit(1.0, "mm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0)) +
  guides(fill = guide_legend(nrow = 1))

p_strip <- (p_left_legged | p_right_legged) +
  plot_layout(widths = c(0.85, 1.0), guides = "keep")

# Width chosen so each subplot column fits its legend without clipping.
W_mm <- 128; H_mm <- 50
agg_png(OUT_PNG, width = W_mm, height = H_mm, units = "mm", res = 1200,
        background = "white")
print(p_strip)
dev.off()

# AI-editable vector exports (RGB, Arial, no textLength so Illustrator can reflow text).
w_in <- W_mm / 25.4; h_in <- H_mm / 25.4
svglite(OUT_SVG, width = w_in, height = h_in)
print(p_strip)
dev.off()
svg_txt <- readLines(OUT_SVG, warn = FALSE)
svg_txt <- gsub('"Liberation Sans"', '"Arial"', svg_txt, fixed = TRUE)
svg_txt <- gsub(" textLength='[^']*'", "", svg_txt)
svg_txt <- gsub(" lengthAdjust='[^']*'", "", svg_txt)
writeLines(svg_txt, OUT_SVG)

cairo_pdf(OUT_PDF, width = w_in, height = h_in, family = BASE_FONT)
print(p_strip)
dev.off()

img <- readPNG(OUT_PNG)
message("Wrote ", OUT_PNG,
        " | dims: ", ncol(img), "x", nrow(img),
        " | aspect: ", round(ncol(img) / nrow(img), 3))
message("Wrote ", OUT_SVG)
message("Wrote ", OUT_PDF)
