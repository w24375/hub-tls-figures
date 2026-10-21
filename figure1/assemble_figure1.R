#!/usr/bin/env Rscript
# Figure 1 (7-panel) assembly for the Hub-TLS manuscript.
# Layout: a = spatial plate (2 samples x [TLS stage, CCL19, CXCL13]);
#         b = radial gradient; c = TLS-ring chemokine heatmap; d = hub-centric gradient meta forest;
#         e = proximal vs distal dumbbell; f = hub cell-type enrichment; g = chemokine producers.
# Backend: R (ggplot2 + patchwork). Panels b-g are native vector; panel a tiles are raster.
# Source data: figure1_assembled/source_data/. fig1c values transcribed from the project
# ring-heatmap figure; fig1d values transcribed from the spatial evidence report (evidence 1 table).

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(grid); library(png)
  library(dplyr); library(tidyr); library(readr); library(svglite); library(ragg)
})

BASE_FONT <- "Arial"; NATURE_SIZE <- 7
FIG_DIR <- "data/hub_tls_paper_outline/figure1_assembled"
SRC <- file.path(FIG_DIR, "source_data")
TILE_DIR <- file.path(FIG_DIR, "panel_a_tiles"); dir.create(TILE_DIR, showWarnings = FALSE)
NAT <- "data/hub_tls_paper_outline/figure1_assembled/panel_a_tiles_svg"  # legend-enlarged tiles

theme_nature <- function(...) {
  theme_classic(base_size = NATURE_SIZE, base_family = BASE_FONT) +
    theme(axis.line = element_line(linewidth = 0.35, colour = "black"),
          axis.ticks = element_line(linewidth = 0.35, colour = "black"),
          axis.text = element_text(size = NATURE_SIZE, colour = "black"),
          axis.title = element_text(size = NATURE_SIZE, colour = "black"),
          plot.title = element_text(size = NATURE_SIZE, face = "bold", hjust = 0),
          legend.title = element_text(size = NATURE_SIZE - 0.5),
          legend.text = element_text(size = NATURE_SIZE - 0.8),
          legend.key.size = unit(3, "mm"), legend.margin = margin(0, 0, 0, 0),
          panel.grid = element_blank(), plot.margin = margin(2, 2, 2, 2, "mm"), ...)
}

## ---------- 1. Panel a tiles (native vector) ----------
marker_label <- c(tls_stage = "TLS stage", CCL19 = "CCL19", CXCL13 = "CXCL13")

blend_white <- function(hex, a) {
  m <- col2rgb(hex) / 255
  out <- a * m + (1 - a)
  rgb(out[1, ], out[2, ], out[3, ])
}

parse_tile <- function(svg_path) {
  ln  <- readLines(svg_path, warn = FALSE)
  cir <- ln[grepl("<circle", ln)]
  cx   <- as.numeric(sub(".*cx='([0-9.]+)'.*", "\\1", cir))
  cy   <- as.numeric(sub(".*cy='([0-9.]+)'.*", "\\1", cir))
  rr   <- as.numeric(sub(".*r='([0-9.]+)'.*",  "\\1", cir))
  fill <- sub(".*fill: (#[0-9A-Fa-f]{6}).*", "\\1", cir)
  fop  <- suppressWarnings(as.numeric(sub(".*fill-opacity: ([0-9.]+).*", "\\1", cir)))
  fop[is.na(fop)] <- 1
  lin <- ln[grepl("<line", ln)]
  seg <- NULL
  if (length(lin)) {
    seg <- data.frame(
      x1 = as.numeric(sub(".*x1='([0-9.-]+)'.*", "\\1", lin)),
      y1 = as.numeric(sub(".*y1='([0-9.-]+)'.*", "\\1", lin)),
      x2 = as.numeric(sub(".*x2='([0-9.-]+)'.*", "\\1", lin)),
      y2 = as.numeric(sub(".*y2='([0-9.-]+)'.*", "\\1", lin)))
  }
  list(spots = data.frame(cx, cy, r = rr, fill = fill, fop = fop), lines = seg)
}

UNITS_TO_MM <- 0.22  # ~56 mm cell / ~255 svg units

vector_tile <- function(sample, marker) {
  svg_path <- file.path(NAT, paste0(sample, "_", marker, "_HE.svg"))
  ln  <- readLines(svg_path, warn = FALSE)
  tl  <- parse_tile(svg_path)
  sp  <- tl$spots
  leg <- sp[sp$cx > 204, ]          # legend keys (tls_stage only)
  sp  <- sp[sp$cx <= 204, ]         # tissue spots
  sp$yp <- -sp$cy                    # flip y for ggplot

  # Original canvas size from the outermost rect (width/height attrs)
  rect_line <- ln[grepl("<rect", ln) & grepl("width='2", ln)][1]
  W <- as.numeric(sub(".*width='([0-9.]+)'.*",  "\\1", rect_line))
  H <- as.numeric(sub(".*height='([0-9.]+)'.*", "\\1", rect_line))
  # Pad shorter side to square (mimics previous white padding);
  # extend right side so absolute-size legend text is not clipped
  side <- max(W, H)
  x_extra <- if (marker == "tls_stage") 30 else 22
  xlim <- c(-8, side + x_extra)
  ylim <- c(-side, 0)

  p <- ggplot(sp, aes(cx, yp)) +
    geom_point(colour = blend_white(sp$fill, sp$fop),
               size = sp$r * 2 * UNITS_TO_MM, shape = 16, stroke = 0)

  # scale bar (native segments + label)
  if (!is.null(tl$lines) && nrow(tl$lines) > 0) {
    sg <- tl$lines; sg$y1p <- -sg$y1; sg$y2p <- -sg$y2
    p <- p + geom_segment(data = sg, aes(x = x1, y = y1p, xend = x2, yend = y2p),
                          colour = "black", linewidth = 0.5, inherit.aes = FALSE)
    hl <- sg[which.max(abs(sg$x2 - sg$x1)), ]
    p <- p + annotate("text", x = (hl$x1 + hl$x2) / 2, y = hl$y1p + 2, label = "500 µm",
                      size = 6.5 / .pt, family = BASE_FONT, hjust = 0.5, vjust = 0)
  }

  if (marker == "tls_stage" && nrow(leg) > 0) {
    # Place legend in the canvas right margin (original zone x>205)
    leg <- leg[order(leg$cy), ]
    leg$label <- c("E-TLS", "PFL-TLS", "SFL-TLS")
    leg_col <- blend_white(leg$fill, leg$fop)
    lx <- 207; ly <- -95            # title anchor
    item_y <- -(c(108, 124, 140))   # item rows (more vertical spacing)
    p <- p +
      annotate("text", x = lx, y = ly, label = "TLS stage", fontface = "bold",
               size = 7 / .pt, family = BASE_FONT, hjust = 0) +
      annotate("point", x = lx + 2.5, y = item_y, colour = leg_col, size = 1.8) +
      annotate("text", x = lx + 11, y = item_y, label = leg$label,
               size = 6.5 / .pt, family = BASE_FONT, hjust = 0)
  }

  if (marker != "tls_stage") {
    pal <- read_csv(file.path(FIG_DIR, sprintf("colorbar_%s_palette.csv", marker)),
                    show_col_types = FALSE)
    # Colorbar in right margin: x from 207, vertically centered on canvas
    bx <- 207; bw <- 7
    ytop <- -80; ybot <- -150
    n <- nrow(pal)
    ys <- seq(ytop, ybot, length.out = n + 1)
    cb <- data.frame(ymax = ys[-(n + 1)], ymin = ys[-1], fill = pal$hex)
    p <- p +
      geom_rect(data = cb, aes(xmin = bx, xmax = bx + bw, ymin = ymin, ymax = ymax,
                               fill = fill), inherit.aes = FALSE) +
      scale_fill_identity() +
      annotate("text", x = bx + bw / 2, y = ytop + 14, label = marker_label[[marker]],
               fontface = "bold", size = 7 / .pt, family = BASE_FONT) +
      annotate("text", x = bx + bw + 3, y = seq(ytop, ybot, length.out = 5),
               label = 4:0, size = 6.5 / .pt, family = BASE_FONT, hjust = 0)
  }

  p + coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE) +
    theme_void() + theme(plot.margin = margin(3, 0.5, 0.5, 0.5, "mm"))
}

text_carrier <- function(txt, rot = 0, size = 7.5, bold = TRUE) {
  ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = txt, angle = rot,
             size = size / .pt, family = BASE_FONT, fontface = if (bold) "bold" else "plain") +
    theme(plot.margin = margin(0, 0, 0, 0, "mm"))
}

hdr <- (text_carrier("") | text_carrier("TLS stage") | text_carrier("CCL19") | text_carrier("CXCL13")) +
  plot_layout(widths = c(0.055, 1, 1, 1))
mk_row <- function(sample, label) {
  (text_carrier(label, rot = 90, size = 7) |
     vector_tile(sample, "tls_stage") |
     vector_tile(sample, "CCL19") |
     vector_tile(sample, "CXCL13")) +
    plot_layout(widths = c(0.055, 1, 1, 1))
}
a_grid <- hdr / mk_row("PDAC_PA05_2", "PDAC PA05-2") / mk_row("Breast_M11", "Breast M11") +
  plot_layout(heights = c(0.07, 1, 1))

## ---------- 2. Panel b: radial gradient ----------
radial_agg <- read_csv(file.path(SRC, "fig1_radial_gradient_source_data.csv"), show_col_types = FALSE) |>
  arrange(dist_bin)
p_b <- ggplot(radial_agg, aes(x = dist_bin, y = mean_hub)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#33B5A5", alpha = 0.25) +
  geom_line(linewidth = 0.8, colour = "#6C8EBF") +
  geom_point(size = 1.5, colour = "#6C8EBF") +
  scale_x_continuous(name = "Distance to nearest TLS (px)", expand = c(0, 0), breaks = seq(0, 1500, 300)) +
  scale_y_continuous(name = "Hub combined score", expand = expansion(mult = c(0.05, 0.1))) +
  annotate("text", x = 1480, y = max(radial_agg$upper) * 0.92,
           label = "Hub signal peaks at TLS\nand decays outward",
           size = NATURE_SIZE * 0.35, hjust = 1, family = BASE_FONT, colour = "#6C8EBF") +
  labs(tag = "B") +
  theme_nature()

## ---------- 3. Panel c: TLS-ring chemokine heatmap ----------
rc <- read_csv(file.path(SRC, "fig1c_tls_ring_chemokine.csv"), show_col_types = FALSE)
ring_levels <- c("TLS-0-100", "TLS-100-200", "TLS-200-300", "TLS-300-500", "TLS-500-1000", "TLS->1000")
rc <- rc |> mutate(ring = factor(ring, levels = ring_levels),
                   gene = factor(gene, levels = c("CCL19", "CCL21", "CXCL10", "CXCL13", "LTB", "TNFSF14")))
p_c <- ggplot(rc, aes(x = ring, y = gene)) +
  geom_tile(aes(fill = mean_expr), colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", mean_expr)), size = NATURE_SIZE * 0.3, family = BASE_FONT) +
  scale_fill_gradient2(low = "#6C8EBF", mid = "white", high = "#C41E3A", midpoint = 0.18,
                       name = "Mean\nlog-norm") +
  labs(x = "Distance from TLS (px)", y = NULL, tag = "C") +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = NATURE_SIZE - 0.8),
        axis.text.y = element_text(size = NATURE_SIZE - 0.7),
        axis.title.x = element_text(size = NATURE_SIZE - 0.5),
        plot.title = element_text(size = NATURE_SIZE - 0.3, face = "bold"),
        legend.title = element_text(size = NATURE_SIZE - 1.2),
        legend.text = element_text(size = NATURE_SIZE - 1.4),
        legend.key.size = unit(2.6, "mm"))

## ---------- 4. Panel d: hub-centric gradient meta forest ----------
meta <- read_csv(file.path(SRC, "fig1d_hub_centric_gradient_meta.csv"), show_col_types = FALSE) |>
  mutate(signature = recode(signature,
                            "CCR7+ T-cell signature" = "CCR7+ T cell",
                            "CXCR5+ B-cell signature" = "CXCR5+ B cell",
                            "Pooled chemokines" = "Pooled chemokines"),
         signature = reorder(signature, -pooled_r),
         plab = ifelse(p_value < 0.001, sprintf("p=%.0e", p_value),
                       ifelse(p_value < 0.01, sprintf("p=%.3f", p_value), sprintf("p=%.2f (ns)", p_value))))
p_d <- ggplot(meta, aes(y = signature, x = pooled_r)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey50") +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.3, linewidth = 0.4, colour = "grey40") +
  geom_point(aes(size = n_samples), colour = "#6C8EBF") +
  geom_text(aes(x = 0.014, label = plab), hjust = 0, size = NATURE_SIZE * 0.3, family = BASE_FONT, colour = "grey25") +
  scale_size_continuous(range = c(1.6, 3.2), name = "Samples") +
  scale_x_continuous(limits = c(-0.095, 0.075), breaks = c(-0.08, -0.04, 0)) +
  labs(x = "Pooled Spearman r (95% CI)", y = NULL, tag = "D") +
  coord_cartesian(clip = "off") +
  theme_nature() +
  theme(axis.text.x = element_text(size = NATURE_SIZE - 1),
        axis.text.y = element_text(size = NATURE_SIZE - 0.8),
        plot.title = element_text(size = NATURE_SIZE - 0.4, face = "bold"),
        plot.margin = margin(2, 10, 2, 2, "mm"),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-4, 0, 0, 0))

## ---------- 5. Panel e: proximal vs distal dumbbell ----------
df <- read_csv(file.path(SRC, "C_proximal_vs_distal_hub_DE_sample_level.csv"), show_col_types = FALSE) |>
  mutate(gene = factor(gene, levels = rev(c("CCL19", "CXCL13", "CCL21", "LTB"))),
         fdr_label = sprintf("FDR = %.0e", fdr),
         pct_label = sprintf("%.0f%%", frac_prox_higher * 100))
long <- df |>
  select(gene, median_prox, median_dist) |>
  pivot_longer(c(median_prox, median_dist), names_to = "group", values_to = "expr") |>
  mutate(group = recode(group, median_dist = "Distal hub (>500 px)", median_prox = "Proximal hub (<200 px)"))
xmax <- max(long$expr)
p_e <- ggplot(long, aes(y = gene)) +
  geom_line(aes(x = expr, group = gene), colour = "grey70", linewidth = 0.5) +
  geom_point(aes(x = expr, colour = group), size = 2.2) +
  geom_text(data = df, aes(x = xmax * 1.12,
            label = sprintf("+%.2f  |  %.0f%%  |  FDR %.0e", median_diff_log2scale, frac_prox_higher * 100, fdr)),
            hjust = 0, size = NATURE_SIZE * 0.34, family = BASE_FONT, colour = "grey20") +
  annotate("text", x = xmax * 1.12, y = 4.62, label = "Diff | % samples | paired FDR",
           hjust = 0, size = NATURE_SIZE * 0.3, family = BASE_FONT, fontface = "bold") +
  scale_colour_manual(values = c("Distal hub (>500 px)" = "grey60", "Proximal hub (<200 px)" = "#C41E3A"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  coord_cartesian(clip = "off") +
  labs(x = "Median log-normalized expression (proximal vs distal, n = 105)", y = NULL, tag = "E") +
  theme_nature() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(6, 30, 2, 2, "mm"))

## ---------- 6. Panel f: hub cell-type enrichment (top hits) ----------
enr <- read_csv(file.path(SRC, "fig1f_hub_celltype_enrichment.csv"), show_col_types = FALSE)
top <- bind_rows(
  enr |> filter(direction == "enriched") |> slice_max(mean_log2FC, n = 6, with_ties = FALSE),
  enr |> filter(direction == "depleted") |> slice_min(mean_log2FC, n = 4, with_ties = FALSE)
) |> mutate(display_name = recode(display_name,
                                  "Fibroblasts (CCL19/CCL21+ FRC)" = "FRC+ fibroblasts",
                                  "Tfh-like / Exhausted T cells" = "Tfh-like/Exh T",
                                  "CD8 T cells (cytotoxic effector)" = "CD8 T (cytotoxic)",
                                  "NK cells (cytotoxic)" = "NK (cytotoxic)",
                                  "Activated T cells (ICOS+)" = "Activated T (ICOS+)",
                                  "Normal epithelial cells" = "Normal epi.",
                                  "Malignant epithelial cells" = "Malignant epi."),
            display_name = reorder(display_name, mean_log2FC))
p_f <- ggplot(top, aes(y = display_name, x = mean_log2FC, fill = direction)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = stars, x = mean_log2FC + ifelse(mean_log2FC > 0, 0.09, -0.09),
                hjust = ifelse(mean_log2FC > 0, 0, 1)),
            size = NATURE_SIZE * 0.36, family = BASE_FONT, colour = "grey25") +
  scale_fill_manual(values = c(enriched = "#1B9E77", depleted = "#D95F02")) +
  scale_x_continuous(limits = c(-1.45, 2.15)) +
  annotate("text", x = -0.72, y = 7.5, label = "Enriched\nin hubs", colour = "#1B9E77",
           size = NATURE_SIZE * 0.30, family = BASE_FONT, fontface = "bold", hjust = 0.5, vjust = 0.5, lineheight = 0.95) +
  annotate("text", x = 1.07, y = 2.5, label = "Depleted\nin hubs", colour = "#D95F02",
           size = NATURE_SIZE * 0.30, family = BASE_FONT, fontface = "bold", hjust = 0.5, vjust = 0.5, lineheight = 0.95) +
  labs(x = "Mean log2 fold-change (hub / non-hub)", y = NULL, tag = "F") +
  theme_nature() +
  theme(axis.text.y = element_text(size = NATURE_SIZE - 0.7), legend.position = "none")

## ---------- 7. Panel g: chemokine producers ----------
prod <- read_csv(file.path(SRC, "fig1g_chemokine_producers.csv"), show_col_types = FALSE)
short <- c("DCs (mregDC LAMP3/CCR7/IDO1+)" = "mregDC",
           "NK cells (cytotoxic, IFN-resp)" = "NK (cytotoxic)",
           "Fibroblasts (CCL19/CCL21+ FRC)" = "FRC+ fibroblast",
           "DCs" = "DCs", "Fibroblasts" = "Fibroblast",
           "Macrophages (CXCL10/CXCL9+ M1-like)" = "M1-like macro.",
           "Endothelial cells" = "Endothelial", "Epithelial cells" = "Epithelial",
           "Malignant epithelial cells (basal/squamous)" = "Malignant epi.",
           "Tfh-like / Exhausted T cells" = "Tfh-like/Exh T",
           "Low-quality / Dying cells" = "Low-quality",
           "Activated/Stressed T cells" = "Activated T",
           "CD8 T cells (IFN-gamma+ effector)" = "CD8 T (IFNg+)")
prod <- prod |> mutate(cell_short = short[cell_type],
                       gene = factor(gene, levels = c("CCL19", "CCL21", "CXCL13")))
ord <- prod |> group_by(cell_short) |> summarise(m = max(avg_counts_log1p)) |> arrange(m) |> pull(cell_short)
p_g <- ggplot(prod, aes(x = gene, y = factor(cell_short, levels = rev(ord)))) +
  geom_point(aes(size = pct_expressing, fill = avg_counts_log1p), shape = 21, colour = "grey40", stroke = 0.3) +
  scale_fill_gradient(low = "white", high = "#2171B5", name = "log1p(avg\ncounts)") +
  scale_size_continuous(range = c(0.8, 3.6), name = "%\nexpressing") +
  labs(x = NULL, y = NULL, tag = "G") +
  theme_nature() +
  theme(axis.text.y = element_text(size = NATURE_SIZE - 1.2),
        axis.text.x = element_text(size = NATURE_SIZE - 1, angle = 30, hjust = 1),
        legend.position = "bottom", legend.direction = "vertical",
        legend.box = "horizontal", legend.box.margin = margin(-3, 0, 0, 0))

## ---------- 8. Compose ----------
a_tag <- ggplot() + theme_void() + labs(tag = "A") +
  theme(plot.tag = element_text(size = 13, face = "bold", family = BASE_FONT),
        plot.tag.position = c(0, 1), plot.margin = margin(0, 0, 0, 0, "mm"))
a_block <- (a_tag | a_grid) + plot_layout(widths = c(0.012, 1))
row2 <- p_b + p_c + p_d + plot_layout(widths = c(1.12, 1.32, 1.0))
row3 <- p_e + p_f + p_g + plot_layout(widths = c(1.42, 1.0, 1.18))
final <- a_block / row2 / row3 + plot_layout(heights = c(2.5, 1, 1.1)) &
  theme(plot.tag = element_text(size = 13, face = "bold", family = BASE_FONT))

W <- 183; H <- 225
w_in <- W / 25.4; h_in <- H / 25.4
svg_out <- file.path(FIG_DIR, "fig1_composite.svg")
svglite(svg_out, width = w_in, height = h_in); print(final); dev.off()
# svglite writes the resolved system family ("Liberation Sans"); relabel as Arial
# (Liberation Sans is metrically identical to Arial; renders as true Arial on machines that have it)
svg_txt <- readLines(svg_out, warn = FALSE)
svg_txt <- gsub('"Liberation Sans"', '"Arial"', svg_txt, fixed = TRUE)
# strip AI-incompatible attrs so reopened SVG stays editable
svg_txt <- gsub(" textLength='[^']*'", "", svg_txt)
svg_txt <- gsub(" lengthAdjust='[^']*'", "", svg_txt)
writeLines(svg_txt, svg_out)
cairo_pdf(file.path(FIG_DIR, "fig1_composite.pdf"), width = w_in, height = h_in, family = BASE_FONT); print(final); dev.off()
ragg::agg_png(file.path(FIG_DIR, "fig1_composite.png"), width = W, height = H, units = "mm", res = 400); print(final); dev.off()
message("Figure 1 written: ", W, " x ", H, " mm")

## ---------- 9. Export individual panels for the assembly report ----------
REP_DIR <- file.path(FIG_DIR, "report_panels"); dir.create(REP_DIR, showWarnings = FALSE)
save_panel <- function(plot, name, w_mm, h_mm) {
  ragg::agg_png(file.path(REP_DIR, name), width = w_mm, height = h_mm, units = "mm", res = 300)
  print(plot); dev.off()
}
save_panel(a_grid, "panel_a_spatial_plate.png", 175, 118)
save_panel(p_b, "panel_b_radial_gradient.png", 90, 55)
save_panel(p_c, "panel_c_ring_heatmap.png", 80, 60)
save_panel(p_d, "panel_d_meta_forest.png", 75, 60)
save_panel(p_e, "panel_e_proximal_distal.png", 95, 60)
save_panel(p_f, "panel_f_celltype_enrichment.png", 80, 62)
save_panel(p_g, "panel_g_chemokine_producers.png", 85, 62)
message("Report panels exported to ", REP_DIR)
