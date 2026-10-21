#!/usr/bin/env Rscript
# Figure 5 (7-panel) assembly for the Hub-TLS manuscript.
# Claim: CCR7+ T / CXCR5+ B cells statistically mediate hub-to-TLS formation;
#        hub proximity licenses the E-to-PFL organizational transition (spatial necessity).
# a = mediation DAG schematic; b = mediation proportions; c = L-R spatial correlation;
# d = L-R spatial proximity (observed/expected); e = proximal vs distal TLS PFL module (n=36);
# f = TLS stage composition by hub proximity; g = module profiles by hub proximity (spot-level).
# Backend: R (ggplot2 + patchwork). All panels native vector.

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(grid)
  library(dplyr); library(tidyr); library(readr); library(svglite); library(ragg)
})

BASE_FONT <- "Arial"; NATURE_SIZE <- 7
FIG_DIR <- "data/hub_tls_paper_outline/figure5_assembled"
SRC <- file.path(FIG_DIR, "source_data")
COL_T <- "#D24B40"; COL_B <- "#66B2A4"; COL_PROX <- "#C41E3A"; COL_DIST <- "grey60"

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

## ---------- 1. Panel a: mediation DAG schematic ----------
box <- function(x0, x1, y0, y1, fill, border, label, size = NATURE_SIZE * 0.36, bold = TRUE) {
  list(annotate("rect", xmin = x0, xmax = x1, ymin = y0, ymax = y1,
                fill = fill, colour = border, linewidth = 0.5),
       annotate("text", x = (x0 + x1) / 2, y = (y0 + y1) / 2, label = label,
                size = size, family = BASE_FONT, fontface = if (bold) "bold" else "plain",
                lineheight = 1.0))
}
arr <- function(x0, y0, x1, y1) {
  annotate("segment", x = x0, y = y0, xend = x1, yend = y1,
           arrow = arrow(length = unit(1.6, "mm"), type = "closed"),
           colour = "grey25", linewidth = 0.5)
}
p_a <- ggplot() +
  coord_cartesian(xlim = c(0, 11), ylim = c(-0.9, 6.2), clip = "off") +
  box(0.1, 3.7, 1.8, 3.8, "#EAF2F8", "#6C8EBF",
      "Hub chemokine signal\n(CCL19/CCL21/CXCL13)", size = NATURE_SIZE * 0.28) +
  box(4.7, 7.3, 3.7, 5.1, "#FDF0EE", COL_T, "CCR7+ T-cell\nrecruitment", size = NATURE_SIZE * 0.27) +
  box(4.7, 7.3, 0.5, 1.9, "#EFF7F5", COL_B, "CXCR5+ B-cell\nrecruitment", size = NATURE_SIZE * 0.27) +
  box(8.4, 10.95, 1.8, 3.8, "#F5F5F5", "grey40",
      "TLS presence\n& maturation", size = NATURE_SIZE * 0.27) +
  arr(3.65, 3.25, 4.75, 4.15) + arr(3.65, 2.35, 4.75, 1.45) +
  arr(7.25, 4.25, 8.38, 3.35) + arr(7.25, 1.35, 8.38, 2.25) +
  annotate("text", x = 2.85, y = 4.55, label = "CCL19/21 -> CCR7",
           size = NATURE_SIZE * 0.3, family = BASE_FONT, colour = "grey30") +
  annotate("text", x = 2.85, y = 0.9, label = "CXCL13 -> CXCR5",
           size = NATURE_SIZE * 0.3, family = BASE_FONT, colour = "grey30") +
  annotate("text", x = 9.7, y = 5.75, hjust = 0.5, label = "10.8% of TLS\nmaturation (p=0.04)",
           size = NATURE_SIZE * 0.3, family = BASE_FONT, colour = COL_T, lineheight = 0.95) +
  annotate("text", x = 9.7, y = 0.15, hjust = 0.5, label = "34.4% of TLS\npresence (p=0.01)",
           size = NATURE_SIZE * 0.3, family = BASE_FONT, colour = "#3D8A80", lineheight = 0.95) +
  annotate("text", x = 5.5, y = -0.75,
           label = "Cross-sectional statistical mediation; causal direction requires functional validation",
           size = NATURE_SIZE * 0.3, family = BASE_FONT, colour = "grey45", fontface = "italic") +
  labs(tag = "A") + theme_void() +
  theme(plot.tag = element_text(size = 13, face = "bold", family = BASE_FONT),
        plot.margin = margin(2, 2, 2, 2, "mm"))

## ---------- 2. Panel b: mediation proportions ----------
med <- read_csv(file.path(SRC, "fig5b_mediation_summary.csv"), show_col_types = FALSE) |>
  mutate(prop_pct = prop_mediated * 100,
         sig = case_when(ACME_p < 0.05 ~ "*", TRUE ~ "ns"),
         mediator = factor(mediator, levels = c("CXCR5+ B", "CCR7+ T")),
         outcome = factor(outcome, levels = c("TLS", "TLS stage"),
                          labels = c("TLS presence", "TLS maturation stage")),
         over = prop_pct > 100)
p_b <- ggplot(med, aes(x = mediator, y = prop_pct, fill = outcome)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65,
           aes(group = outcome, linetype = over, linewidth = over), colour = "grey20") +
  geom_text(aes(label = sprintf("%.1f%%\n%s", prop_pct, sig), group = outcome),
            position = position_dodge(width = 0.75), vjust = -0.15,
            size = NATURE_SIZE * 0.3, family = BASE_FONT, lineheight = 0.9) +
  annotate("text", x = 0.62, y = 90, label = "over-mediation\n(model unstable)",
           size = NATURE_SIZE * 0.3, family = BASE_FONT, colour = "grey40",
           fontface = "italic", lineheight = 0.9, hjust = 0) +
  annotate("segment", x = 1.30, y = 90, xend = 1.66, yend = 104,
           colour = "grey40", linewidth = 0.3,
           arrow = arrow(length = unit(1, "mm"), type = "closed")) +
  scale_fill_manual(values = c("TLS presence" = "#E08214", "TLS maturation stage" = "#F1C48B"), name = NULL) +
  scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "dashed"), guide = "none") +
  scale_linewidth_manual(values = c(`FALSE` = 0.3, `TRUE` = 0.7), guide = "none") +
  scale_y_continuous(limits = c(0, 145), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Proportion mediated (%)", tag = "B") +
  theme_nature() +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.background = element_blank())

## ---------- 3. Panel c: L-R spatial correlation (full cohort, sample-level) ----------
lr_c <- read_csv(file.path(SRC, "fig5c_lr_correlation_full_cohort.csv"), show_col_types = FALSE) |>
  filter(!is.na(rho_ccr7) & !is.na(rho_cxcr5))
n_c <- nrow(lr_c)
p_pair_c <- wilcox.test(lr_c$rho_ccr7, lr_c$rho_cxcr5, paired = TRUE)$p.value
sig_c <- c(mean(lr_c$p_ccr7 < 0.05 & lr_c$rho_ccr7 > 0),
           mean(lr_c$p_cxcr5 < 0.05 & lr_c$rho_cxcr5 > 0))
lr_long <- lr_c |>
  pivot_longer(c(rho_ccr7, rho_cxcr5), names_to = "pair", values_to = "rho") |>
  mutate(pair = factor(pair, levels = c("rho_ccr7", "rho_cxcr5"),
                       labels = c("CCL19/CCL21\n– CCR7", "CXCL13\n– CXCR5")))
y0_c <- floor(min(lr_long$rho) * 10) / 10 - 0.09
y1_c <- ceiling(max(lr_long$rho) * 10) / 10 + 0.07
p_c <- ggplot(lr_long, aes(x = pair, y = rho)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey55") +
  geom_boxplot(aes(fill = pair), width = 0.52, outlier.shape = NA,
               colour = "grey20", linewidth = 0.35, alpha = 0.9) +
  geom_jitter(aes(colour = pair), width = 0.13, size = 0.6, alpha = 0.45, stroke = 0) +
  annotate("text", x = 1:2, y = y0_c + 0.035, size = NATURE_SIZE * 0.3, family = BASE_FONT,
           label = sprintf("%.0f%% p<0.05", 100 * sig_c), colour = "grey30") +
  annotate("text", x = 1.5, y = y1_c, vjust = 1, size = NATURE_SIZE * 0.3, family = BASE_FONT,
           label = sprintf("n = %d samples\npaired Wilcoxon p = %.1e", n_c, p_pair_c), lineheight = 0.95) +
  scale_fill_manual(values = c("CCL19/CCL21\n– CCR7" = "#6C8EBF", "CXCL13\n– CXCR5" = "#66B2A4")) +
  scale_colour_manual(values = c("CCL19/CCL21\n– CCR7" = "#6C8EBF", "CXCL13\n– CXCR5" = "#66B2A4")) +
  scale_y_continuous(limits = c(y0_c, y1_c + 0.02)) +
  labs(x = NULL, y = "Spearman rho (hub spots)", tag = "C") +
  theme_nature() +
  theme(legend.position = "none", axis.text.x = element_text(size = NATURE_SIZE - 1, lineheight = 0.9))

## ---------- 4. Panel d: L-R spatial proximity (full cohort, sample-level) ----------
lrp <- read_csv(file.path(SRC, "fig5d_lr_proximity_full_cohort.csv"), show_col_types = FALSE) |>
  mutate(pair = factor(pair, levels = c("CCL19/CCL21 – CCR7", "CXCL13 – CXCR5"),
                       labels = c("CCL19/CCL21\n– CCR7", "CXCL13\n– CXCR5")))
d_stats <- lrp |> group_by(pair) |>
  summarise(n = n(), pct_sig = 100 * mean(perm_p < 0.05),
            p = wilcox.test(log(ratio))$p.value, .groups = "drop") |>
  mutate(stars = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns"))
y1_d <- max(1.25, max(lrp$ratio) * 1.08)
p_d <- ggplot(lrp, aes(x = pair, y = ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4, colour = "grey40") +
  geom_boxplot(aes(fill = pair), width = 0.52, outlier.shape = NA,
               colour = "grey20", linewidth = 0.35, alpha = 0.9) +
  geom_jitter(aes(colour = pair), width = 0.13, size = 0.6, alpha = 0.45, stroke = 0) +
  geom_text(data = d_stats, aes(x = pair, y = y1_d * 0.965, label = stars),
            size = NATURE_SIZE * 0.4, family = BASE_FONT, vjust = 1) +
  geom_text(data = d_stats, aes(x = pair, y = 0.045,
                                label = sprintf("n=%d, %.0f%% perm p<0.05", n, pct_sig)),
            size = NATURE_SIZE * 0.27, family = BASE_FONT, colour = "grey30") +
  scale_fill_manual(values = c("CCL19/CCL21\n– CCR7" = "#6C8EBF", "CXCL13\n– CXCR5" = "#66B2A4")) +
  scale_colour_manual(values = c("CCL19/CCL21\n– CCR7" = "#6C8EBF", "CXCL13\n– CXCR5" = "#66B2A4")) +
  scale_y_continuous(limits = c(0, y1_d)) +
  labs(x = NULL, y = "Observed / expected distance", tag = "D") +
  theme_nature() +
  theme(legend.position = "none", axis.text.x = element_text(size = NATURE_SIZE - 1, lineheight = 0.9))

## ---------- 5. Panel e: proximal vs distal TLS PFL module (sample-level) ----------
pd <- read_csv(file.path(SRC, "fig5e_pfl_tls_proximal_vs_distal_sample_level.csv"), show_col_types = FALSE)
n_e <- nrow(pd); med_diff <- median(pd$diff); frac_pos <- mean(pd$diff > 0)
p_e_val <- wilcox.test(pd$prox, pd$dist, paired = TRUE)$p.value
lim_e <- c(0, 0.58)
p_e <- ggplot(pd, aes(x = dist, y = prox)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey50") +
  geom_point(size = 1.4, alpha = 0.8, colour = COL_PROX) +
  annotate("text", x = 0.572, y = 0.02, hjust = 1, vjust = 0, family = BASE_FONT, size = NATURE_SIZE * 0.32,
           label = sprintf("n = %d samples\nmedian diff = +%.3f\npaired Wilcoxon p = %.1e\n%.0f%% samples prox > dist",
                           n_e, med_diff, p_e_val, frac_pos * 100), lineheight = 0.95) +
  scale_x_continuous(limits = lim_e) +
  scale_y_continuous(limits = lim_e) +
  labs(x = "PFL-TLS module, hub-distal TLS\n(≥500 px ≈ 230 µm)",
       y = "PFL-TLS module, hub-proximal TLS\n(≤200 px ≈ 90 µm)", tag = "E") +
  theme_nature() +
  theme(aspect.ratio = 1,
        axis.title = element_text(size = NATURE_SIZE - 0.5, lineheight = 0.95))

## ---------- 6. Panel f: TLS stage composition by hub proximity (sample-level) ----------
d_spots <- read_csv(file.path(SRC, "D_tls_spot_dist_to_hub.csv"), show_col_types = FALSE)
f_stats <- read_csv(file.path(SRC, "fig5f_stage_composition_sample_level.csv"), show_col_types = FALSE)
star_of <- function(f) case_when(f < 0.001 ~ "***", f < 0.01 ~ "**", f < 0.05 ~ "*", TRUE ~ "ns")
star_map <- setNames(star_of(f_stats$fdr), f_stats$stage)

per_s <- d_spots |>
  count(sample, grp, tls_stage) |>
  pivot_wider(names_from = tls_stage, values_from = n, values_fill = 0) |>
  mutate(tot = `E-TLS` + `PFL-TLS` + `SFL-TLS`,
         across(c(`E-TLS`, `PFL-TLS`, `SFL-TLS`), ~ .x / tot)) |>
  select(-tot)
keep_f <- per_s |>
  pivot_wider(names_from = grp, values_from = c(`E-TLS`, `PFL-TLS`, `SFL-TLS`), names_glue = "{grp}_{.value}") |>
  # group retained iff it contributed >=5 spots; reconstruct spot counts
  inner_join(
    d_spots |> count(sample, grp, name = "tot") |>
      pivot_wider(names_from = grp, values_from = tot),
    by = "sample") |>
  filter(proximal >= 5 & distal >= 5) |> pull(sample)
n_f <- length(keep_f)
comp_f <- per_s |>
  filter(sample %in% keep_f) |>
  pivot_longer(c(`E-TLS`, `PFL-TLS`, `SFL-TLS`), names_to = "tls_stage", values_to = "frac") |>
  group_by(grp, tls_stage) |>
  summarise(frac = mean(frac), .groups = "drop") |>
  mutate(pct = 100 * frac,
         hub_proximity = factor(grp, levels = c("distal", "proximal"),
                                labels = c("Hub-distal TLS", "Hub-proximal TLS")),
         tls_stage = factor(tls_stage, levels = c("SFL-TLS", "PFL-TLS", "E-TLS")),
         star = star_map[as.character(tls_stage)],
         seg_lab = sprintf("%.1f%%%s", pct, ifelse(star == "ns", " ns", paste0(" ", star))))
p_f <- ggplot(comp_f, aes(x = hub_proximity, y = pct, fill = tls_stage)) +
  geom_col(width = 0.62, colour = "grey20", linewidth = 0.3) +
  geom_text(aes(label = seg_lab), position = position_stack(vjust = 0.5),
            size = NATURE_SIZE * 0.28, family = BASE_FONT) +
  annotate("text", x = 1.5, y = 106.5, size = NATURE_SIZE * 0.27, family = BASE_FONT,
           colour = "grey30", lineheight = 0.95,
           label = sprintf("n = %d paired samples; mean of per-sample fractions\nstars: paired Wilcoxon FDR", n_f)) +
  scale_fill_manual(values = c("E-TLS" = "#D24B40", "PFL-TLS" = "#E69F00", "SFL-TLS" = "#5DA05D"), name = "TLS stage") +
  scale_x_discrete(labels = c("Hub-distal TLS" = "Hub-distal TLS\n≥500 px ≈ 230 µm",
                              "Hub-proximal TLS" = "Hub-proximal TLS\n≤200 px ≈ 90 µm")) +
  scale_y_continuous(limits = c(0, 112), breaks = c(0, 25, 50, 75, 100),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Mean % of TLS spots", tag = "F") +
  theme_nature() +
  theme(legend.position = "right", axis.text.x = element_text(size = NATURE_SIZE - 1, lineheight = 0.9))

## ---------- 7. Panel g: module profiles by hub proximity (spot-level) ----------
mod <- read_csv(file.path(SRC, "fig5g_module_by_hub_proximity.csv"), show_col_types = FALSE)
mod_lab <- c("E_TLS_GENES" = "E-TLS genes", "PFL_TLS_GENES" = "PFL-TLS genes",
             "SFL_TLS_GENES" = "SFL-TLS genes", "Hub_chemo" = "Hub chemo",
             "Hub_tcell" = "Hub T-cell", "Hub_stem" = "Hub stem")
mod <- mod |>
  mutate(fdr_num = suppressWarnings(as.numeric(sub("^<", "", as.character(fdr)))),
         fdr_stars = case_when(fdr_num < 0.001 ~ "***", fdr_num < 0.01 ~ "**",
                               fdr_num < 0.05 ~ "*", TRUE ~ "ns"),
         label = factor(mod_lab[module], levels = rev(mod_lab)))
mod_long <- mod |> select(label, median_proximal, median_distal) |>
  pivot_longer(c(median_proximal, median_distal), names_to = "grp", values_to = "score") |>
  mutate(grp = recode(grp, median_proximal = "Hub-proximal TLS", median_distal = "Hub-distal TLS"))
p_g <- ggplot(mod_long, aes(y = label)) +
  geom_hline(yintercept = 3.5, linetype = "dotted", linewidth = 0.3, colour = "grey85") +
  geom_line(aes(x = score, group = label), colour = "grey70", linewidth = 0.5) +
  geom_point(aes(x = score, colour = grp), size = 2) +
  geom_text(data = mod,
            aes(x = max(mod_long$score) * 1.12,
                label = sprintf("%+.3f %s", median_proximal - median_distal, fdr_stars)),
            hjust = 0, size = NATURE_SIZE * 0.3, family = BASE_FONT, colour = "grey20") +
  annotate("text", x = max(mod_long$score) * 1.12, y = 6.7, label = "Prox - dist",
           hjust = 0, size = NATURE_SIZE * 0.3, family = BASE_FONT, fontface = "bold") +
  scale_colour_manual(values = c("Hub-proximal TLS" = COL_PROX, "Hub-distal TLS" = COL_DIST), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  coord_cartesian(clip = "off") +
  labs(x = "Median module score", y = NULL, tag = "G") +
  theme_nature() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(4, 18, 2, 2, "mm"))

## ---------- 8. Compose ----------
row1 <- p_a + p_b + plot_layout(widths = c(1.25, 1))
row2 <- p_c + p_d + p_e + plot_layout(widths = c(1.05, 1.05, 1))
row3 <- p_f + p_g + plot_layout(widths = c(0.95, 1.3))
final <- row1 / row2 / row3 + plot_layout(heights = c(1, 1, 1)) &
  theme(plot.tag = element_text(size = 13, face = "bold", family = BASE_FONT))

W <- 183; H <- 205
w_in <- W / 25.4; h_in <- H / 25.4
svg_out <- file.path(FIG_DIR, "fig5_composite.svg")
svglite(svg_out, width = w_in, height = h_in); print(final); dev.off()
# svglite writes the resolved system family ("Liberation Sans"); relabel as Arial
svg_txt <- readLines(svg_out, warn = FALSE)
svg_txt <- gsub('"Liberation Sans"', '"Arial"', svg_txt, fixed = TRUE)
# strip AI-incompatible attrs so reopened SVG stays editable
svg_txt <- gsub(" textLength='[^']*'", "", svg_txt)
svg_txt <- gsub(" lengthAdjust='[^']*'", "", svg_txt)
writeLines(svg_txt, svg_out)
cairo_pdf(file.path(FIG_DIR, "fig5_composite.pdf"), width = w_in, height = h_in, family = BASE_FONT); print(final); dev.off()
ragg::agg_png(file.path(FIG_DIR, "fig5_composite.png"), width = W, height = H, units = "mm", res = 400); print(final); dev.off()
message("Figure 5 written: ", W, " x ", H, " mm")

## ---------- 9. Export individual panels ----------
REP_DIR <- file.path(FIG_DIR, "report_panels"); dir.create(REP_DIR, showWarnings = FALSE)
save_panel <- function(plot, name, w_mm, h_mm) {
  ragg::agg_png(file.path(REP_DIR, name), width = w_mm, height = h_mm, units = "mm", res = 300)
  print(plot); dev.off()
}
save_panel(p_a, "panel_a_mediation_dag.png", 100, 55)
save_panel(p_b, "panel_b_mediation_proportions.png", 80, 62)
save_panel(p_c, "panel_c_lr_correlation.png", 80, 58)
save_panel(p_d, "panel_d_lr_proximity.png", 80, 58)
save_panel(p_e, "panel_e_proximal_distal_scatter.png", 80, 62)
save_panel(p_f, "panel_f_stage_composition.png", 75, 62)
save_panel(p_g, "panel_g_module_dumbbell.png", 90, 62)
message("Report panels exported.")
