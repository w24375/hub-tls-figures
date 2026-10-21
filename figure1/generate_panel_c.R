#!/usr/bin/env Rscript
# Fig 1 panel c (sample-level): proximal vs distal hub differential expression.
# Data: pseudoreplication_fix sample-level paired Wilcoxon (n = 105 samples).
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr); library(svglite); library(ragg)
})

BASE_FONT <- "Arial"; NATURE_SIZE <- 7
theme_nature <- function(...) {
  theme_classic(base_size = NATURE_SIZE, base_family = BASE_FONT) +
    theme(axis.line = element_line(linewidth = 0.35, colour = "black"),
          axis.ticks = element_line(linewidth = 0.35, colour = "black"),
          axis.text = element_text(size = NATURE_SIZE, colour = "black"),
          axis.title = element_text(size = NATURE_SIZE, colour = "black"),
          plot.title = element_text(size = NATURE_SIZE + 1, face = "bold", hjust = 0),
          legend.title = element_text(size = NATURE_SIZE - 0.5),
          legend.text = element_text(size = NATURE_SIZE - 0.8),
          legend.key.size = unit(3, "mm"), legend.margin = margin(0, 0, 0, 0),
          panel.grid = element_blank(), plot.margin = margin(2, 2, 2, 2, "mm"), ...)
}

OUT_DIR <- "data/hub_tls_paper_outline/figure1_assembled"
df <- read_csv(file.path(OUT_DIR, "source_data", "C_proximal_vs_distal_hub_DE_sample_level.csv"),
               show_col_types = FALSE)

df <- df |>
  mutate(gene = factor(gene, levels = rev(c("CCL19", "CXCL13", "CCL21", "LTB"))),
         fdr_label = sprintf("FDR = %.1e", fdr),
         pct_label = sprintf("%.0f%%", frac_prox_higher * 100))

long <- df |>
  select(gene, median_prox, median_dist) |>
  pivot_longer(cols = c(median_prox, median_dist), names_to = "group", values_to = "expr") |>
  mutate(group = recode(group, median_dist = "Distal hub (>500 px)", median_prox = "Proximal hub (<200 px)"))

p <- ggplot(long, aes(y = gene)) +
  geom_line(aes(x = expr, group = gene), colour = "grey70", linewidth = 0.5) +
  geom_point(aes(x = expr, colour = group), size = 2.2) +
  geom_text(data = df, aes(x = max(long$expr) * 1.12, label = paste0("+", round(median_diff_log2scale, 2))),
            hjust = 0, size = NATURE_SIZE * 0.36, family = BASE_FONT, colour = "#B2182B") +
  geom_text(data = df, aes(x = max(long$expr) * 1.45, label = pct_label),
            hjust = 0, size = NATURE_SIZE * 0.36, family = BASE_FONT, colour = "grey30") +
  geom_text(data = df, aes(x = max(long$expr) * 1.75, label = fdr_label),
            hjust = 0, size = NATURE_SIZE * 0.36, family = BASE_FONT, colour = "grey30") +
  annotate("text", x = max(long$expr) * c(1.12, 1.45, 1.75), y = 4.62,
           label = c("Median\ndiff", "Samples\nprox>dist", "Paired\nFDR"),
           hjust = 0, size = NATURE_SIZE * 0.32, family = BASE_FONT, fontface = "bold", lineheight = 0.9) +
  scale_colour_manual(values = c("Distal hub (>500 px)" = "grey60", "Proximal hub (<200 px)" = "#B2182B"),
                      name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  coord_cartesian(clip = "off") +
  labs(title = "Proximal vs distal hubs (sample-level, n = 105)",
       x = "Median log-normalized expression (per-sample means)", y = NULL) +
  theme_nature() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(8, 32, 2, 2, "mm"))

w_in <- 120 / 25.4; h_in <- 70 / 25.4
svglite(file.path(OUT_DIR, "panel_c_proximal_vs_distal_sample_level.svg"), width = w_in, height = h_in); print(p); dev.off()
cairo_pdf(file.path(OUT_DIR, "panel_c_proximal_vs_distal_sample_level.pdf"), width = w_in, height = h_in, family = BASE_FONT); print(p); dev.off()
ragg::agg_png(file.path(OUT_DIR, "panel_c_proximal_vs_distal_sample_level.png"), width = 120, height = 70, units = "mm", res = 600); print(p); dev.off()
message("done")
