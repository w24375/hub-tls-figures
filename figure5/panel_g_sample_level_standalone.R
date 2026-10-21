#!/usr/bin/env Rscript
# Standalone sample-level Panel G for Fig 5 (AI-editable output).
# Per-sample median module scores in hub-proximal vs hub-distal TLS (n = 36 paired
# samples); right column: median difference + paired-Wilcoxon FDR stars.
# Data: fig5g_module_by_hub_proximity_sample_level.csv (sample-level, canonical modules).

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr); library(svglite); library(ragg)
})

BASE_FONT <- "Arial"; NATURE_SIZE <- 7
SRC <- "data/hub_tls_paper_outline/figure5_assembled/source_data"
OUT <- "data/hub_tls_paper_outline/figures_AI_ready"
COL_PROX <- "#C41E3A"; COL_DIST <- "grey60"

theme_nature <- function(...) {
  theme_classic(base_size = NATURE_SIZE, base_family = BASE_FONT) +
    theme(axis.line = element_line(linewidth = 0.35, colour = "black"),
          axis.ticks = element_line(linewidth = 0.35, colour = "black"),
          axis.text = element_text(size = NATURE_SIZE, colour = "black"),
          axis.title = element_text(size = NATURE_SIZE, colour = "black"),
          legend.title = element_text(size = NATURE_SIZE - 0.5),
          legend.text = element_text(size = NATURE_SIZE - 0.8),
          legend.key.size = unit(3, "mm"), legend.margin = margin(0, 0, 0, 0),
          panel.grid = element_blank(), plot.margin = margin(4, 20, 2, 2, "mm"), ...)
}

mod <- read_csv(file.path(SRC, "fig5g_module_by_hub_proximity_sample_level.csv"),
                show_col_types = FALSE)
lab_map <- c("E_TLS_GENES" = "E-TLS genes", "PFL_TLS_GENES" = "PFL-TLS genes",
             "SFL_TLS_GENES" = "SFL-TLS genes", "Hub_chemo" = "Hub chemo",
             "Hub T-cell" = "Hub T-cell", "Hub_stem" = "Hub stem")
star_of <- function(f) ifelse(f < 0.001, "***", ifelse(f < 0.01, "**", ifelse(f < 0.05, "*", "ns")))
mod <- mod |>
  mutate(label = factor(lab_map[module], levels = rev(unname(lab_map))),
         stars = star_of(fdr),
         diff_lab = sprintf("%+.3f %s", median_diff, stars))

mod_long <- mod |>
  select(label, median_prox, median_dist) |>
  pivot_longer(c(median_prox, median_dist), names_to = "grp", values_to = "score") |>
  mutate(grp = recode(grp, median_prox = "Hub-proximal TLS", median_dist = "Hub-distal TLS"))

x_lab <- max(mod_long$score) * 1.10
p_g <- ggplot(mod_long, aes(y = label)) +
  geom_hline(yintercept = 3.5, linetype = "dotted", linewidth = 0.3, colour = "grey85") +
  geom_line(aes(x = score, group = label), colour = "grey70", linewidth = 0.5) +
  geom_point(aes(x = score, colour = grp), size = 2) +
  geom_text(data = mod, aes(x = x_lab, label = diff_lab),
            hjust = 0, size = NATURE_SIZE * 0.3, family = BASE_FONT, colour = "grey20") +
  annotate("text", x = x_lab, y = 6.65, label = "Prox - dist", hjust = 0,
           size = NATURE_SIZE * 0.3, family = BASE_FONT, fontface = "bold") +
  annotate("text", x = x_lab, y = 0.35, hjust = 0, size = NATURE_SIZE * 0.26,
           family = BASE_FONT, colour = "grey40", lineheight = 0.95,
           label = "paired Wilcoxon FDR\nn = 36 samples") +
  scale_colour_manual(values = c("Hub-proximal TLS" = COL_PROX, "Hub-distal TLS" = COL_DIST), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  coord_cartesian(clip = "off") +
  labs(x = "Median module score (per-sample medians)", y = NULL, tag = "G") +
  theme_nature() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.tag = element_text(size = 13, face = "bold", family = BASE_FONT))

W <- 90; H <- 66
svg_out <- file.path(OUT, "panel_G_module_dumbbell_sample_level.svg")
svglite(svg_out, width = W / 25.4, height = H / 25.4); print(p_g); dev.off()
svg_txt <- readLines(svg_out, warn = FALSE)
svg_txt <- gsub('"Liberation Sans"', '"Arial"', svg_txt, fixed = TRUE)
svg_txt <- gsub(" textLength='[^']*'", "", svg_txt)
svg_txt <- gsub(" lengthAdjust='[^']*'", "", svg_txt)
writeLines(svg_txt, svg_out)
cairo_pdf(file.path(OUT, "panel_G_module_dumbbell_sample_level.pdf"),
          width = W / 25.4, height = H / 25.4, family = BASE_FONT); print(p_g); dev.off()
ragg::agg_png(file.path(OUT, "panel_G_module_dumbbell_sample_level.png"),
              width = W, height = H, units = "mm", res = 400); print(p_g); dev.off()
message("panel G sample-level written to ", OUT)
