#!/usr/bin/env Rscript
# Figure 2 (7-panel) assembly for the Hub-TLS manuscript.
# Claim: Stem hubs lie closer to TLS than Base hubs (directional maturation axis).
# Layout: a = spatial maps (PDAC + Breast, Base/Stem/TLS);
#         b = per-sample paired scatter (n=138); c = log2FC histogram; d = per-cancer forest;
#         e = proximity-threshold bars (spot-level); f = module dumbbell (n=132);
#         g = cellular composition difference (n=156).
# Backend: R (ggplot2 + patchwork). Panels b-g native vector; panel a raster tiles.
# All sample-level stats from source_data/ (fig2a = updated sample-level pipeline, 2026-07-25).

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(grid); library(png)
  library(dplyr); library(tidyr); library(readr); library(svglite); library(ragg)
})

BASE_FONT <- "Arial"; NATURE_SIZE <- 7
FIG_DIR <- "data/hub_tls_paper_outline/figure2_assembled"
SRC <- file.path(FIG_DIR, "source_data")
TILE_DIR <- file.path(FIG_DIR, "panel_a_tiles"); dir.create(TILE_DIR, showWarnings = FALSE)
NAT <- "data/spatial/hub_precursor_analysis/signal_promotion/figures_nature"
COL_BASE <- "#6C8EBF"; COL_STEM <- "#C41E3A"

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

cancer_type_from_dataset <- function(dataset) {
  case_when(
    grepl("Breast", dataset, ignore.case = TRUE) ~ "Breast",
    grepl("Lung", dataset, ignore.case = TRUE) ~ "Lung",
    grepl("Ovarian", dataset, ignore.case = TRUE) ~ "Ovarian",
    grepl("PDAC|Pancrea", dataset, ignore.case = TRUE) ~ "PDAC",
    grepl("GBM|Glioblastoma", dataset, ignore.case = TRUE) ~ "GBM",
    grepl("HCC|Hepatocellular", dataset, ignore.case = TRUE) ~ "HCC",
    grepl("Gastric|Stomach", dataset, ignore.case = TRUE) ~ "Gastric",
    grepl("Colorectal|Colon|Rectal", dataset, ignore.case = TRUE) ~ "Colorectal",
    grepl("Prostate", dataset, ignore.case = TRUE) ~ "Prostate",
    grepl("Bladder", dataset, ignore.case = TRUE) ~ "Bladder",
    grepl("Renal", dataset, ignore.case = TRUE) ~ "Renal",
    grepl("Oral", dataset, ignore.case = TRUE) ~ "Oral",
    TRUE ~ "Other")
}

## ---------- 1. Panel a tiles: directionality proof maps ----------
process_tile <- function(sample) {
  svg_path <- file.path(NAT, paste0(sample, "_directionality_proof_HE.svg"))
  svg_txt <- paste(readLines(svg_path, warn = FALSE), collapse = "\n")
  svg_txt <- gsub(paste0("<text[^>]*>", sample, "[^<]*</text>"), "", svg_txt)
  mod_svg <- tempfile(fileext = ".svg"); writeLines(svg_txt, mod_svg)
  tmp <- tempfile(fileext = ".png")
  system2("rsvg-convert", c("-w", "1200", shQuote(mod_svg), "-o", shQuote(tmp)), stdout = FALSE, stderr = FALSE)
  arr <- png::readPNG(tmp)
  if (length(dim(arr)) == 2) arr <- array(rep(arr, 3), c(dim(arr), 3))
  if (dim(arr)[3] == 4) { a <- arr[, , 4]; for (ch in 1:3) arr[, , ch] <- arr[, , ch] * a + (1 - a); arr <- arr[, , 1:3] }
  gray <- apply(arr, c(1, 2), mean)
  ink <- gray < 0.97
  rows_ink <- which(rowMeans(ink) > 0.0005); cols_ink <- which(colMeans(ink) > 0.0005)
  r0 <- max(1, min(rows_ink) - 12); r1 <- min(nrow(arr), max(rows_ink) + 12)
  c0 <- max(1, min(cols_ink) - 12); c1 <- min(ncol(arr), max(cols_ink) + 12)
  arr <- arr[r0:r1, c0:c1, , drop = FALSE]
  h2 <- nrow(arr); w2 <- ncol(arr); s <- max(h2, w2)
  new <- array(1, c(s, s, 3))
  ro <- floor((s - h2) / 2); co <- floor((s - w2) / 2)
  new[(ro + 1):(ro + h2), (co + 1):(co + w2), ] <- arr
  out <- file.path(TILE_DIR, paste0(sample, "__directionality.png"))
  png::writePNG(new, out)
  out
}
samples_main <- c("PDAC_PA05_2", "Breast_M11")
tile_files <- sapply(samples_main, process_tile)

raster_carrier <- function(png_path) {
  img <- png::readPNG(png_path)
  g <- rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc"))
  ggplot() + annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    coord_cartesian(clip = "off") + theme_void() +
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"), aspect.ratio = 1)
}
text_carrier <- function(txt, rot = 0, size = 7.5, bold = TRUE) {
  ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = txt, angle = rot,
             size = size / .pt, family = BASE_FONT, fontface = if (bold) "bold" else "plain") +
    theme(plot.margin = margin(0, 0, 0, 0, "mm"))
}

## ---------- 2. Load data ----------
dist_s <- read_csv(file.path(SRC, "fig2a_source.csv"), show_col_types = FALSE) |>
  mutate(cancer = cancer_type_from_dataset(dataset))
modules <- read_csv(file.path(SRC, "fig2b_stem_vs_base_sample_level.csv"), show_col_types = FALSE)
prox <- read_csv(file.path(SRC, "hub_subtype_distance_to_tls_summary.csv"), show_col_types = FALSE)
comp <- read_csv(file.path(SRC, "fig2g_base_stem_composition_full_cohort.csv"), show_col_types = FALSE)

n_s <- nrow(dist_s); n_closer <- sum(dist_s$direction == "Stem closer")
pct_closer <- round(100 * n_closer / n_s, 1)
wp <- wilcox.test(dist_s$stem_mean, dist_s$base_mean, paired = TRUE)$p.value
mean_lfc <- mean(dist_s$log2FC); ci_l <- mean_lfc - 1.96 * sd(dist_s$log2FC) / sqrt(n_s)
ci_u <- mean_lfc + 1.96 * sd(dist_s$log2FC) / sqrt(n_s)

## ---------- 3. Panel b: paired scatter ----------
p_b <- ggplot(dist_s, aes(x = base_mean, y = stem_mean, colour = cancer)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey50") +
  geom_point(size = 1.1, alpha = 0.85) +
  annotate("text", x = 30, y = 2600, hjust = 0, vjust = 1, family = BASE_FONT, size = NATURE_SIZE * 0.34,
           label = sprintf("%d/%d (%.1f%%) samples:\nStem hub closer\npaired Wilcoxon p = %.1e",
                           n_closer, n_s, pct_closer, wp), lineheight = 0.95) +
  scale_x_log10(name = "Base hub mean distance to TLS (px)") +
  scale_y_log10(name = "Stem hub mean distance to TLS (px)") +
  scale_colour_brewer(palette = "Paired", name = "Cancer") +
  labs(tag = "B") +
  theme_nature() +
  theme(legend.position = "right", legend.key.size = unit(2.4, "mm"),
        legend.text = element_text(size = NATURE_SIZE - 1.2))

## ---------- 4. Panel c: log2FC histogram ----------
dist_s <- dist_s |>
  mutate(sig_cat = factor(signif_cat,
                          levels = c("Stem closer (FDR<0.05)", "Not significant", "Base closer (FDR<0.05)")))
p_c <- ggplot(dist_s, aes(x = log2FC, fill = sig_cat)) +
  geom_histogram(binwidth = 0.12, colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey40") +
  annotate("text", x = -2.2, y = 15.8, hjust = 0, vjust = 1, family = BASE_FONT, size = NATURE_SIZE * 0.34,
           label = sprintf("mean = %.2f [%.2f, %.2f]\nStouffer Z = -48.98", mean_lfc, ci_l, ci_u),
           lineheight = 0.95) +
  scale_fill_manual(values = c("Stem closer (FDR<0.05)" = "#33B5A5",
                               "Not significant" = "grey70",
                               "Base closer (FDR<0.05)" = "#D95F02"), name = NULL) +
  labs(x = "log2FC (Stem / Base mean distance)", y = "Number of samples", tag = "C") +
  theme_nature() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.text = element_text(size = NATURE_SIZE - 1.4),
        legend.key.size = unit(2.2, "mm")) +
  guides(fill = guide_legend(nrow = 1))

## ---------- 5. Panel d: per-cancer forest (dot strip + median) ----------
per_cancer <- dist_s |>
  group_by(cancer) |>
  summarise(n = n(), med = median(log2FC), .groups = "drop") |>
  arrange(med) |>
  mutate(cancer = factor(cancer, levels = cancer))
dist_s2 <- dist_s |> mutate(cancer = factor(cancer, levels = levels(per_cancer$cancer)))
p_d <- ggplot(dist_s2, aes(y = cancer, x = log2FC)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey40") +
  geom_point(position = position_jitter(height = 0.18, seed = 1), size = 0.7, alpha = 0.55, colour = "#6C8EBF") +
  geom_point(data = per_cancer, aes(x = med, y = cancer), shape = 18, size = 2.6, colour = COL_STEM) +
  geom_text(data = per_cancer, aes(x = 1.62, label = paste0("n=", n)), hjust = 0,
            size = NATURE_SIZE * 0.32, family = BASE_FONT, colour = "grey30") +
  scale_x_continuous(limits = c(-2.6, 2.2)) +
  labs(x = "log2FC (Stem / Base mean distance)", y = NULL, tag = "D") +
  theme_nature() +
  theme(plot.title = element_text(size = NATURE_SIZE - 0.3, face = "bold"))

## ---------- 6. Panel e: proximity thresholds (spot-level, descriptive) ----------
prox_l <- prox |>
  select(hub_subtype, n, pct_within_100, pct_within_200, pct_within_500) |>
  pivot_longer(starts_with("pct_within"), names_to = "threshold", values_to = "pct") |>
  mutate(threshold = factor(sub("pct_within_", "<=", threshold),
                            levels = c("<=100", "<=200", "<=500")),
         hub_subtype = factor(hub_subtype, levels = c("Base hub", "Stem hub")))
p_e <- ggplot(prox_l, aes(x = threshold, y = pct, fill = hub_subtype)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  geom_text(aes(label = sprintf("%.1f%%", pct), group = hub_subtype),
            position = position_dodge(width = 0.75), vjust = -0.4,
            size = NATURE_SIZE * 0.3, family = BASE_FONT) +
  scale_fill_manual(values = c("Base hub" = COL_BASE, "Stem hub" = COL_STEM), name = NULL) +
  scale_y_continuous(limits = c(0, 92), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Distance threshold to nearest TLS (px)", y = "% of hub spots", tag = "E") +
  theme_nature() +
  theme(legend.position = c(0.28, 0.8), legend.background = element_blank(),
        axis.text.x = element_text(size = NATURE_SIZE - 0.6),
        plot.title = element_text(size = NATURE_SIZE - 0.4, face = "bold"))

## ---------- 7. Panel f: module dumbbell (Stem vs Base, n=132) ----------
mod_lab <- c("chemokine_score" = "Chemokine", "myeloid_score" = "Myeloid",
             "tls_maturity_score" = "TLS maturity", "Hub_tcell" = "Hub T-cell",
             "Hub_dc" = "Hub DC", "HEV_GENES" = "HEV genes",
             "E_TLS_GENES" = "E-TLS genes", "PFL_TLS_GENES" = "PFL-TLS genes",
             "SFL_TLS_GENES" = "SFL-TLS genes")
mf <- modules |>
  mutate(label = mod_lab[module],
         note = sprintf("+%.3f | %.0f%% | FDR %.0e", median_diff, frac_stem_higher * 100, fdr))
mf$note[mf$module == "myeloid_score"] <-
  sprintf("+%.3f | %.0f%% | FDR %.1e (weak)", mf$median_diff[mf$module == "myeloid_score"],
          mf$frac_stem_higher[mf$module == "myeloid_score"] * 100, mf$fdr[mf$module == "myeloid_score"])
mf <- mf |> mutate(label = reorder(label, median_diff))
mf_long <- mf |> select(label, median_base, median_stem) |>
  pivot_longer(c(median_base, median_stem), names_to = "subtype", values_to = "score") |>
  mutate(subtype = recode(subtype, median_base = "Base hub", median_stem = "Stem hub"))
xmax_f <- max(mf_long$score)
p_f <- ggplot(mf_long, aes(y = label)) +
  geom_line(aes(x = score, group = label), colour = "grey70", linewidth = 0.5) +
  geom_point(aes(x = score, colour = subtype), size = 2) +
  geom_text(data = mf, aes(x = xmax_f * 1.10, label = note),
            hjust = 0, size = NATURE_SIZE * 0.3, family = BASE_FONT, colour = "grey20") +
  annotate("text", x = xmax_f * 1.10, y = 9.7, label = "Diff | % samples | paired FDR",
           hjust = 0, size = NATURE_SIZE * 0.3, family = BASE_FONT, fontface = "bold") +
  scale_colour_manual(values = c("Base hub" = COL_BASE, "Stem hub" = COL_STEM), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  coord_cartesian(clip = "off") +
  labs(x = "Median module score (per-sample means)", y = NULL, tag = "F") +
  theme_nature() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(4, 34, 2, 2, "mm"))

## ---------- 8. Panel g: composition difference ----------
comp2 <- comp |>
  mutate(stem_enrichment = -mean_log2FC,
         stars = case_when(fdr < 0.001 ~ "***", fdr < 0.01 ~ "**", fdr < 0.05 ~ "*", TRUE ~ ""),
         display = recode(cell_type,
                          "NK cells (cytotoxic, IFN-resp)" = "NK cells (cytotoxic)",
                          "Tfh-like / Exhausted T cells" = "Tfh-like / Exhausted T",
                          "Activated T cells (ICOS+)" = "Activated T (ICOS+)",
                          "CD4 T cells (naive/TCF7+)" = "CD4 T (naive/TCF7+)",
                          "DCs (mregDC LAMP3/CCR7/IDO1+)" = "mregDC",
                          "CD8 T cells (cytotoxic effector)" = "CD8 T (cytotoxic)",
                          "Low-quality / Dying cells" = "Low-quality / Dying",
                          "DCs (mature DC)" = "DCs (mature)",
                          "Macrophages (CXCL10/CXCL9+ M1-like)" = "M1-like macrophages",
                          "Malignant epithelial cells (basal/squamous)" = "Malignant epi. (basal)"))
top_g <- bind_rows(
  comp2 |> slice_max(stem_enrichment, n = 6, with_ties = FALSE),
  comp2 |> slice_min(stem_enrichment, n = 4, with_ties = FALSE)
) |> mutate(display = reorder(display, stem_enrichment),
            dir = ifelse(stem_enrichment > 0, "Higher in Stem", "Higher in Base"))
p_g <- ggplot(top_g, aes(y = display, x = stem_enrichment, fill = dir)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = stars, x = stem_enrichment + ifelse(stem_enrichment > 0, 0.12, -0.12),
                hjust = ifelse(stem_enrichment > 0, 0, 1)),
            size = NATURE_SIZE * 0.36, family = BASE_FONT, colour = "grey25") +
  scale_fill_manual(values = c("Higher in Stem" = COL_STEM, "Higher in Base" = COL_BASE)) +
  scale_x_continuous(limits = c(min(top_g$stem_enrichment) - 0.7, max(top_g$stem_enrichment) + 0.8)) +
  labs(x = "Mean log2 fold-change (Stem / Base)", y = NULL, tag = "G") +
  theme_nature() +
  theme(axis.text.y = element_text(size = NATURE_SIZE - 0.5), legend.position = "none")

## ---------- 9. Compose ----------
a_tag <- ggplot() + theme_void() + labs(tag = "A") +
  theme(plot.tag = element_text(size = 13, face = "bold", family = BASE_FONT),
        plot.tag.position = c(0, 1), plot.margin = margin(0, 0, 0, 0, "mm"))
a_maps <- (text_carrier("PDAC PA05-2") | text_carrier("Breast M11")) /
  (raster_carrier(tile_files["PDAC_PA05_2"]) | raster_carrier(tile_files["Breast_M11"])) +
  plot_layout(heights = c(0.06, 1))
a_row <- (a_tag | a_maps | p_b) + plot_layout(widths = c(0.02, 2, 1.25))
row2 <- p_c + p_d + p_e + plot_layout(widths = c(1.1, 1, 0.95))
row3 <- p_f + p_g + plot_layout(widths = c(1.4, 1))
final <- a_row / row2 / row3 + plot_layout(heights = c(1, 0.95, 1.05)) &
  theme(plot.tag = element_text(size = 13, face = "bold", family = BASE_FONT))

W <- 183; H <- 195
w_in <- W / 25.4; h_in <- H / 25.4
svg_out <- file.path(FIG_DIR, "fig2_composite.svg")
svglite(svg_out, width = w_in, height = h_in); print(final); dev.off()
# svglite writes the resolved system family ("Liberation Sans"); relabel as Arial
svg_txt <- readLines(svg_out, warn = FALSE)
svg_txt <- gsub('"Liberation Sans"', '"Arial"', svg_txt, fixed = TRUE)
# strip AI-incompatible attrs so reopened SVG stays editable
svg_txt <- gsub(" textLength='[^']*'", "", svg_txt)
svg_txt <- gsub(" lengthAdjust='[^']*'", "", svg_txt)
writeLines(svg_txt, svg_out)
cairo_pdf(file.path(FIG_DIR, "fig2_composite.pdf"), width = w_in, height = h_in, family = BASE_FONT); print(final); dev.off()
ragg::agg_png(file.path(FIG_DIR, "fig2_composite.png"), width = W, height = H, units = "mm", res = 400); print(final); dev.off()
message("Figure 2 written: ", W, " x ", H, " mm; samples: ", n_s, ", stem closer: ", n_closer, " (", pct_closer, "%)")

## ---------- 10. Export individual panels for the report ----------
REP_DIR <- file.path(FIG_DIR, "report_panels"); dir.create(REP_DIR, showWarnings = FALSE)
save_panel <- function(plot, name, w_mm, h_mm) {
  ragg::agg_png(file.path(REP_DIR, name), width = w_mm, height = h_mm, units = "mm", res = 300)
  print(plot); dev.off()
}
save_panel((text_carrier("PDAC PA05-2") | text_carrier("Breast M11")) /
             (raster_carrier(tile_files["PDAC_PA05_2"]) | raster_carrier(tile_files["Breast_M11"])) +
             plot_layout(heights = c(0.06, 1)), "panel_a_spatial_maps.png", 120, 62)
save_panel(p_b, "panel_b_paired_scatter.png", 90, 62)
save_panel(p_c, "panel_c_log2fc_histogram.png", 80, 58)
save_panel(p_d, "panel_d_per_cancer_forest.png", 75, 58)
save_panel(p_e, "panel_e_proximity_bars.png", 70, 58)
save_panel(p_f, "panel_f_module_dumbbell.png", 100, 62)
save_panel(p_g, "panel_g_composition_diff.png", 80, 62)
message("Report panels exported.")
