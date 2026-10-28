#!/usr/bin/env Rscript
# Paired slope plot: per-sample mean continuous-slingshot pseudotime,
# Hub-only vs TLS-only, on the 131 samples that have all three classes
# (matching the manuscript's 96.9% / p < 2e-16 figure, which used means).
# The Transition class is omitted because its per-sample position is unstable
# (mean PT 9.46 > TLS-only 8.28).

library(ggplot2); library(dplyr); library(tidyr); library(readr)
library(svglite); library(ragg)

source(if (file.exists("code/utils/theme.R")) "code/utils/theme.R" else "utils/theme.R")

SRC <- "data/hub_tls_paper_outline/figure3_assembled/source_data"
OUT <- "data/hub_tls_paper_outline/figure3_assembled/report_panels"

cs <- read_csv(file.path(SRC, "fig3j_continuous_slingshot_spots.csv"), show_col_types = FALSE)
complete3 <- cs |>
  filter(!is.na(pseudotime)) |>
  distinct(orig.ident, continuous_class) |>
  count(orig.ident) |>
  filter(n == 3) |>
  pull(orig.ident)
med <- cs |>
  filter(continuous_class %in% c("Hub_only", "TLS_only"), !is.na(pseudotime),
         orig.ident %in% complete3) |>
  group_by(orig.ident, continuous_class) |>
  summarise(med_pt = mean(pseudotime), .groups = "drop") |>
  pivot_wider(names_from = continuous_class, values_from = med_pt) |>
  drop_na()
stopifnot(nrow(med) == 131)

wil <- wilcox.test(med$Hub_only, med$TLS_only, paired = TRUE)
frac_up <- mean(med$Hub_only < med$TLS_only)
n_up <- sum(med$Hub_only < med$TLS_only)
meds <- med |> summarise(across(c(Hub_only, TLS_only), median))
message(sprintf("Hub-first: %d/%d (%.1f%%), paired Wilcoxon p = %.2e",
                n_up, nrow(med), 100 * frac_up, wil$p.value))

long <- med |>
  mutate(up = Hub_only < TLS_only) |>
  pivot_longer(c(Hub_only, TLS_only), names_to = "cls", values_to = "pt") |>
  mutate(x = as.numeric(factor(cls, levels = c("Hub_only", "TLS_only"))))
med_line <- data.frame(x = 1:2, pt = unlist(meds))

p <- ggplot(long, aes(x = x, y = pt)) +
  geom_line(aes(group = orig.ident, colour = up), linewidth = 0.35, alpha = 0.35) +
  geom_point(aes(colour = up), size = 0.7, alpha = 0.5) +
  geom_line(data = med_line, aes(x = x, y = pt), colour = "black", linewidth = 1.1) +
  geom_point(data = med_line, aes(x = x, y = pt), size = 2.6, colour = "black") +
  annotate("label", x = 1.5, y = max(long$pt), vjust = 1.1, hjust = 0.5,
           size = 6 / .pt, family = BASE_FONT, colour = "grey25", lineheight = 0.95,
           fill = "white", alpha = 0.9, label.padding = unit(0.15, "lines"),
           label = sprintf("Hub-first in %d of %d samples (%.1f%%)\npaired Wilcoxon p = %.1e",
                           n_up, nrow(med), 100 * frac_up, wil$p.value)) +
  annotate("text", x = c(1, 2), y = unlist(meds), vjust = -0.8, hjust = c(0.9, 0.1),
           size = 6 / .pt, family = BASE_FONT, fontface = "bold",
           label = sprintf("%.2f", unlist(meds))) +
  scale_colour_manual(values = c("TRUE" = "#6C8EBF", "FALSE" = "#C41E3A"), guide = "none") +
  scale_x_continuous(breaks = 1:2, labels = c("Hub-only", "TLS-only"),
                     limits = c(0.75, 2.25)) +
  labs(x = NULL, y = "Mean slingshot pseudotime\n(per sample)") +
  theme_nature()

W <- 62; H <- 72
sv <- file.path(OUT, "panel_slingshot_paired_slope_standalone.svg")
svglite(sv, width = W / 25.4, height = H / 25.4); print(p); dev.off()
txt <- readLines(sv, warn = FALSE)
txt <- gsub('"Liberation Sans"', '"Arial"', txt, fixed = TRUE)
txt <- gsub(" textLength='[^']*'", "", txt); txt <- gsub(" lengthAdjust='[^']*'", "", txt)
writeLines(txt, sv)
cairo_pdf(file.path(OUT, "panel_slingshot_paired_slope_standalone.pdf"),
          width = W / 25.4, height = H / 25.4, family = BASE_FONT); print(p); dev.off()
agg_png(file.path(OUT, "panel_slingshot_paired_slope_standalone.png"),
        width = W, height = H, units = "mm", res = 600); print(p); dev.off()
message("paired slope panel written")
