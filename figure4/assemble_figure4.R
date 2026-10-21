#!/usr/bin/env Rscript
# Figure 4 (10-panel) assembly for the Hub-TLS manuscript.
# Claim: CCL19/CCL21/CXCL13 are the top-3 genes required to maintain hub and TLS
# signals, confirmed by two independent virtual-KO methods (CellOracle, GEARS),
# with leave-one-out de-circularization and combination-KO collapse.
# Panels: a = impact ranking (both methods); b = method consistency (r = 0.944);
#         c = module-delta heatmap (CellOracle); d = LOO de-circularization dumbbell;
#         e = TLS_like state shift per KO; f = combination-KO collapse;
#         g = PCA displacement field (raster); h = spatial perturbation maps (raster);
#         i = LOO state-collapse bars; j = state-dissolution series (raster).
# Backend: R (ggplot2 + patchwork). Panels a-f,i native vector; g/h/j raster embeds.

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(grid); library(png)
  library(dplyr); library(tidyr); library(readr); library(svglite); library(ragg)
})

BASE_FONT <- "Arial"; NATURE_SIZE <- 7
FIG_DIR <- "data/hub_tls_paper_outline/figure4_assembled"
SRC <- file.path(FIG_DIR, "source_data")
TILE_DIR <- file.path(FIG_DIR, "panel_rasters"); dir.create(TILE_DIR, showWarnings = FALSE)
VK_FIG <- "data/spatial/hub_precursor_analysis/virtual_knockout/figures/flashy"

COL_CO <- "#6C8EBF"; COL_GE <- "#E7A03D"; COL_TRIAD <- "#C41E3A"

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

trim_png <- function(in_path, out_path) {
  arr <- png::readPNG(in_path)
  if (length(dim(arr)) == 2) arr <- array(rep(arr, 3), c(dim(arr), 3))
  if (dim(arr)[3] == 4) { a <- arr[, , 4]; for (ch in 1:3) arr[, , ch] <- arr[, , ch] * a + (1 - a); arr <- arr[, , 1:3] }
  gray <- apply(arr, c(1, 2), mean); ink <- gray < 0.97
  rows_ink <- which(rowMeans(ink) > 0.0005); cols_ink <- which(colMeans(ink) > 0.0005)
  arr <- arr[max(1, min(rows_ink) - 8):min(nrow(arr), max(rows_ink) + 8),
             max(1, min(cols_ink) - 8):min(ncol(arr), max(cols_ink) + 8), , drop = FALSE]
  png::writePNG(arr, out_path); out_path
}
raster_panel <- function(png_path, tag = NULL, margin_mm = 1) {
  img <- png::readPNG(png_path); asp <- nrow(img) / ncol(img)
  p <- ggplot() +
    annotation_custom(rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc")),
                      xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    coord_cartesian(clip = "off") + theme_void() +
    theme(plot.margin = margin(margin_mm, margin_mm, margin_mm, margin_mm, "mm"), aspect.ratio = asp)
  if (!is.null(tag)) p <- p + labs(tag = tag)
  p
}

## ---------- 1. Load data ----------
mc <- read_csv(file.path(SRC, "fig4ab_ko_method_comparison.csv"), show_col_types = FALSE)
mod_delta <- read_csv(file.path(SRC, "fig4c_module_delta_celloracle.csv"), show_col_types = FALSE)
loo <- read_csv(file.path(SRC, "fig4d_loo_decircularization.csv"), show_col_types = FALSE)
shift <- read_csv(file.path(SRC, "fig4e_state_shift.csv"), show_col_types = FALSE)
combo <- read_csv(file.path(SRC, "fig4f_combination_ko.csv"), show_col_types = FALSE)
co_combo <- read_csv(file.path(SRC, "fig4f_celloracle_combination_ko.csv"), show_col_types = FALSE)
ge_single <- read_csv(file.path(SRC, "fig4k_gears_single_ko.csv"), show_col_types = FALSE)

TRIAD <- c("CCL19", "CCL21", "CXCL13")
mc <- mc |> mutate(mean_impact = (impact_score_CO + impact_score_GEARS) / 2,
                   is_triad = gene %in% TRIAD)
gene_order <- mc$gene[order(mc$mean_impact, decreasing = TRUE)]
r_cons <- cor(mc$impact_score_CO, mc$impact_score_GEARS)
message("VERIFY method consistency r = ", round(r_cons, 3), " (expected 0.944)")

## ---------- 2. Panel a: impact ranking ----------
mc_long <- mc |>
  select(gene, impact_score_CO, impact_score_GEARS, is_triad) |>
  pivot_longer(c(impact_score_CO, impact_score_GEARS), names_to = "method", values_to = "impact") |>
  mutate(method = recode(method, impact_score_CO = "CellOracle", impact_score_GEARS = "GEARS"),
         gene = factor(gene, levels = rev(gene_order)))
p_a <- ggplot(mc_long, aes(y = gene, x = impact, fill = method)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66) +
  scale_fill_manual(values = c("CellOracle" = COL_CO, "GEARS" = COL_GE), name = NULL) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5), labels = c("0", "0.25", "0.5"),
                     expand = expansion(mult = c(0, 0.06))) +
  labs(x = "Impact score", y = NULL, tag = "A") +
  theme_nature() +
  theme(axis.text.y = element_text(face = ifelse(rev(gene_order) %in% TRIAD, "bold", "plain"),
                                   colour = ifelse(rev(gene_order) %in% TRIAD, COL_TRIAD, "black")),
        axis.text.x = element_text(size = NATURE_SIZE - 1.2),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.box.margin = margin(-4, 0, 0, 0))

## ---------- 3. Panel b: method consistency ----------
p_b <- ggplot(mc, aes(x = impact_score_CO, y = impact_score_GEARS)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey50") +
  geom_point(aes(colour = is_triad, size = is_triad), show.legend = FALSE) +
  ggrepel::geom_text_repel(data = mc |> filter(is_triad),
                           aes(label = gene), colour = COL_TRIAD, size = NATURE_SIZE * 0.34,
                           family = BASE_FONT, show.legend = FALSE, seed = 7,
                           min.segment.length = 0.2, segment.linewidth = 0.25,
                           max.overlaps = 20, point.padding = 0.6) +
  annotate("text", x = 0.62, y = 0.02, hjust = 1, size = NATURE_SIZE * 0.36, family = BASE_FONT,
           label = sprintf("r = %.3f", r_cons), colour = "grey20", lineheight = 0.95) +
  scale_colour_manual(values = c("TRUE" = COL_TRIAD, "FALSE" = "grey45")) +
  scale_size_manual(values = c("TRUE" = 2.4, "FALSE" = 1.4)) +
  scale_x_continuous(expand = expansion(mult = c(0.16, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  labs(x = "CellOracle impact score", y = "GEARS impact score", tag = "B") +
  theme_nature()

## ---------- 4. Panel c: module-delta heatmap (CellOracle) ----------
md <- mod_delta |>
  select(gene, Hub = Hub_delta_mean, TLS = TLS_delta_mean, Chemo_core = Chemo_core_delta_mean,
         `B cell` = Bcell_delta_mean, `T cell` = Tcell_delta_mean, HEV = HEV_delta_mean) |>
  pivot_longer(-gene, names_to = "module", values_to = "delta") |>
  mutate(gene = factor(gene, levels = rev(gene_order)),
         module = factor(module, levels = c("Hub", "TLS", "Chemo_core", "B cell", "T cell", "HEV")))
p_c <- ggplot(md, aes(x = module, y = gene)) +
  geom_tile(aes(fill = delta), colour = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(delta == 0, "0", sprintf("%.2f", delta))),
            size = NATURE_SIZE * 0.26, family = BASE_FONT) +
  scale_fill_gradient2(low = "#C41E3A", mid = "white", high = "#6C8EBF", midpoint = 0,
                       limits = c(-0.42, 0.42), name = "Module\ndelta") +
  labs(x = NULL, y = NULL, tag = "C") +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = ifelse(rev(gene_order) %in% TRIAD, "bold", "plain"),
                                   colour = ifelse(rev(gene_order) %in% TRIAD, COL_TRIAD, "black")))

## ---------- 5. Panel d: LOO de-circularization ----------
loo2 <- loo |>
  mutate(ym = factor(paste0(gene, " · ", module),
                     levels = rev(paste0(rep(TRIAD, each = 3), " · ", rep(c("Chemo_core", "Hub", "TLS"), 3)))),
         pct_retained = round(100 * loo_mean / orig_mean))
loo_long <- loo2 |>
  select(ym, orig_mean, loo_mean) |>
  pivot_longer(c(orig_mean, loo_mean), names_to = "version", values_to = "delta") |>
  mutate(version = recode(version, orig_mean = "Original KO", loo_mean = "After LOO de-circularization"))
p_d <- ggplot(loo_long, aes(y = ym)) +
  geom_line(aes(x = delta, group = ym), colour = "grey70", linewidth = 0.5,
            arrow = arrow(length = unit(1.2, "mm"), type = "closed")) +
  geom_point(aes(x = delta, colour = version), size = 2) +
  geom_text(data = loo2, aes(x = 0.012, label = sprintf("%d%%", pct_retained)),
            hjust = 0, size = NATURE_SIZE * 0.28, family = BASE_FONT, colour = "grey25") +
  annotate("text", x = 0.012, y = 9.85, label = "LOO retained", hjust = 0,
           size = NATURE_SIZE * 0.28, family = BASE_FONT, fontface = "bold", colour = "grey25") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey50") +
  scale_colour_manual(values = c("Original KO" = "grey45", "After LOO de-circularization" = COL_GE),
                      name = NULL) +
  scale_x_continuous(limits = c(-0.40, 0.10)) +
  coord_cartesian(clip = "off") +
  labs(x = "Module delta after KO", y = NULL, tag = "D") +
  theme_nature() +
  theme(axis.text.y = element_text(size = NATURE_SIZE - 1.2),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.text = element_text(size = NATURE_SIZE - 1.2),
        plot.margin = margin(8, 4, 2, 2, "mm"))

## ---------- 6. Panel e: TLS_like state shift ----------
shift2 <- shift |>
  mutate(gene = factor(gene, levels = rev(gene_order)),
         delta_pp = TLS_like_delta * 100)
shift_long <- shift2 |>
  select(gene, TLS_like_orig, TLS_like_ko) |>
  pivot_longer(c(TLS_like_orig, TLS_like_ko), names_to = "state", values_to = "prop") |>
  mutate(state = recode(state, TLS_like_orig = "Before KO", TLS_like_ko = "After KO"))
p_e <- ggplot(shift_long, aes(y = gene)) +
  geom_line(aes(x = prop, group = gene), colour = "grey70", linewidth = 0.5) +
  geom_point(aes(x = prop, colour = state), size = 2) +
  geom_text(data = shift2, aes(x = 0.545, label = sprintf("%+.1f pp", delta_pp)),
            hjust = 0, size = NATURE_SIZE * 0.28, family = BASE_FONT, colour = "grey25") +
  scale_colour_manual(values = c("Before KO" = "grey50", "After KO" = COL_CO), name = NULL) +
  scale_x_continuous(limits = c(0.25, 0.62)) +
  coord_cartesian(clip = "off") +
  labs(x = "TLS-like spot proportion", y = NULL, tag = "E") +
  theme_nature() +
  theme(axis.text.y = element_text(face = ifelse(rev(gene_order) %in% TRIAD, "bold", "plain"),
                                   colour = ifelse(rev(gene_order) %in% TRIAD, COL_TRIAD, "black"),
                                   size = NATURE_SIZE - 0.8),
        legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(2, 12, 2, 2, "mm"))

## ---------- 7. Panel f: single vs combination KO collapse, both methods ----------
co_f <- co_combo |>
  select(perturbation = combination, collapse_pct = TLS_like_to_Other_percent) |>
  mutate(method = "CellOracle")
ge_f <- bind_rows(
  ge_single |> filter(gene %in% c("CCL19", "CCL21")) |>
    transmute(perturbation = gene, collapse_pct = -TLS_like_delta * 100, method = "GEARS"),
  combo |> transmute(perturbation = combination, collapse_pct = TLS_like_to_Other_percent,
                     method = "GEARS"))
both_f <- bind_rows(co_f, ge_f) |>
  filter(perturbation %in% c("CCL19", "CCL21", "CCL19+CCL21", "CCL19+CCL21+CXCL13")) |>
  mutate(perturbation = factor(perturbation,
                               levels = c("CCL19", "CCL21", "CCL19+CCL21", "CCL19+CCL21+CXCL13"),
                               labels = c("CCL19", "CCL21", "CCL19+21", "Triple")),
         method = factor(method, levels = c("CellOracle", "GEARS")))
# additive expectation = sum of the two single KOs, drawn as dashed segments
# centred on the double-KO bars (category 3; dodge 0.72 -> bar centres +/- 0.18)
add_seg <- data.frame(x0 = c(3 - 0.18, 3 + 0.18),
                      y0 = c(20.72 + 11.99, 4.72 + 4.37))
p_f <- ggplot(both_f, aes(x = perturbation, y = collapse_pct, fill = method)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  geom_text(aes(label = sprintf("%.1f", collapse_pct), group = method),
            position = position_dodge(width = 0.72), vjust = -0.4,
            size = NATURE_SIZE * 0.26, family = BASE_FONT, colour = "black") +
  annotate("segment", x = add_seg$x0 - 0.31, xend = add_seg$x0 + 0.31,
           y = add_seg$y0, yend = add_seg$y0,
           linetype = "dashed", linewidth = 0.5, colour = "black") +
  scale_fill_manual(values = c("CellOracle" = COL_CO, "GEARS" = COL_GE), name = NULL) +
  scale_x_discrete(expand = expansion(add = c(0.6, 1.0))) +
  scale_y_continuous(limits = c(0, 50), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "TLS-like \u2192 Other (%)", tag = "F") +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = NATURE_SIZE - 1.8),
        plot.title = element_text(size = NATURE_SIZE - 0.3, face = "bold"),
        legend.position = "bottom", legend.direction = "horizontal") +
  annotate("segment", x = 0.40, xend = 0.72, y = 48.0, yend = 48.0,
           linetype = "dashed", linewidth = 0.5, colour = "black") +
  annotate("text", x = 0.78, y = 48.0, hjust = 0, vjust = 0.5, size = NATURE_SIZE * 0.26,
           family = BASE_FONT, colour = "grey35", lineheight = 0.9,
           label = "= additive")

## ---------- 8. Panel i: LOO state-transition bars (native) ----------
loo_st <- read_csv(file.path(SRC, "fig4j_loo_state_transition.csv"), show_col_types = FALSE) |>
  mutate(gene = factor(gene, levels = c("CCL19", "CCL21", "CXCL13")),
         method = recode(method, orig = "Original KO", loo = "After LOO de-circularization"),
         pct = tls_to_other * 100)
p_i <- ggplot(loo_st, aes(x = gene, y = pct, fill = method)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66) +
  geom_text(aes(label = sprintf("%.1f%%", pct), group = method),
            position = position_dodge(width = 0.72), vjust = -0.4,
            size = NATURE_SIZE * 0.32, family = BASE_FONT) +
  scale_fill_manual(values = c("Original KO" = "grey55", "After LOO de-circularization" = "#E7A03D"),
                    name = NULL) +
  scale_y_continuous(limits = c(0, 11.5), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "LOO: state collapse",
       x = NULL, y = "TLS-like \u2192 Other (%)", tag = "I") +
  theme_nature() +
  theme(axis.text.x = element_text(size = NATURE_SIZE - 1.0),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.text = element_text(size = NATURE_SIZE - 1.4))

## ---------- 8b. Panel h: spatial perturbation maps (TLS module, LOO, enhanced) ----------
sp35 <- trim_png(file.path(VK_FIG, "ko_spatial_tls_delta_enhanced.png"),
                 file.path(TILE_DIR, "spatial_tls_delta_loo_enhanced.png"))
p_h <- raster_panel(sp35, tag = "H", margin_mm = 2)

## ---------- 8c. Panel g: PCA displacement field (enhanced raster) ----------
pca_png <- trim_png(file.path(VK_FIG, "ko_pca_field_enhanced.png"),
                    file.path(TILE_DIR, "pca_vector_field_enhanced.png"))
p_g <- raster_panel(pca_png, tag = "G", margin_mm = 2)

## ---------- 8d. Panel j: state dissolution series (enhanced raster) ----------
diss_png <- trim_png(file.path(VK_FIG, "ko_state_dissolution_enhanced.png"),
                     file.path(TILE_DIR, "state_dissolution.png"))
p_j <- raster_panel(diss_png, tag = "J", margin_mm = 2)

## ---------- 9. Compose ----------
design <- "
  aaabbbcccccc
  dddddeeeefff
  iiiigggggggg
  hhhhhhhhhhhh
  jjjjjjjjjjjj
"
final <- wrap_plots(p_a, p_b, p_c, p_d, p_e, p_f, p_g, p_h, p_i, p_j, design = design) +
  plot_layout(heights = c(1.0, 1.0, 0.9, 0.9, 0.82)) &
  theme(plot.tag = element_text(size = 13, face = "bold", family = BASE_FONT))

W <- 183; H <- 247
w_in <- W / 25.4; h_in <- H / 25.4
svg_out <- file.path(FIG_DIR, "fig4_composite.svg")
svglite(svg_out, width = w_in, height = h_in); print(final); dev.off()
# svglite writes the resolved system family ("Liberation Sans"); relabel as Arial
svg_txt <- readLines(svg_out, warn = FALSE)
svg_txt <- gsub('"Liberation Sans"', '"Arial"', svg_txt, fixed = TRUE)
# strip AI-incompatible attrs so reopened SVG stays editable
svg_txt <- gsub(" textLength='[^']*'", "", svg_txt)
svg_txt <- gsub(" lengthAdjust='[^']*'", "", svg_txt)
writeLines(svg_txt, svg_out)
cairo_pdf(file.path(FIG_DIR, "fig4_composite.pdf"), width = w_in, height = h_in, family = BASE_FONT); print(final); dev.off()
ragg::agg_png(file.path(FIG_DIR, "fig4_composite.png"), width = W, height = H, units = "mm", res = 400); print(final); dev.off()
message("Figure 4 written: ", W, " x ", H, " mm")

## ---------- 10. Export individual panels ----------
REP_DIR <- file.path(FIG_DIR, "report_panels"); dir.create(REP_DIR, showWarnings = FALSE)
save_panel <- function(plot, name, w_mm, h_mm) {
  ragg::agg_png(file.path(REP_DIR, name), width = w_mm, height = h_mm, units = "mm", res = 300)
  print(plot); dev.off()
}
save_panel(p_a, "panel_a_impact_ranking.png", 75, 62)
save_panel(p_b, "panel_b_method_consistency.png", 62, 62)
save_panel(p_c, "panel_c_module_delta_heatmap.png", 75, 62)
save_panel(p_d, "panel_d_loo_dumbbell.png", 80, 66)
save_panel(p_e, "panel_e_state_shift.png", 75, 62)
save_panel(p_f, "panel_f_combination_ko.png", 70, 50)
save_panel(p_g, "panel_g_pca_field.png", 95, 62)
save_panel(p_h, "panel_h_spatial_perturbation.png", 175, 46)
save_panel(p_i, "panel_i_loo_state_transition.png", 70, 55)
save_panel(p_j, "panel_j_state_dissolution.png", 183, 42)
message("Report panels exported.")
