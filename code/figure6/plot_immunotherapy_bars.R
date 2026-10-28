#!/usr/bin/env Rscript
# Immunotherapy cohort figures (n = 19)
# Fig A: post-treatment TLS status vs clinical response
# Fig B: pre-treatment Hub status vs post-treatment TLS / response
# Fig C: Hub -> TLS -> response chain (alluvial)
# Cross-tab:
#   Hub+ (n=14): TLS+/R 13, TLS-/NR 1
#   Hub- (n=5):  TLS+/R 1, TLS+/NR 2, TLS-/R 1, TLS-/NR 1
library(ggplot2)

outdir <- "data/hub_tls_paper_outline"
fam <- "Arial"

c_red  <- "#A6514A"   # TLS+ (muted brick; replaced #C41E3A as too glaring)
c_blue <- "#6C8EBF"   # response
c_grey <- "#D9D9D9"   # absent / not achieved

theme_set(
  theme_classic(base_size = 6.5, base_family = fam) +
    theme(
      axis.line       = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks      = element_line(linewidth = 0.35, colour = "black"),
      axis.text       = element_text(colour = "black", size = 6.2),
      axis.title.y    = element_text(size = 6.5),
      legend.title    = element_blank(),
      legend.text     = element_text(size = 5.8),
      legend.key.size = unit(3, "mm"),
      legend.position = "top",
      legend.margin   = margin(0, 0, 1, 0),
      panel.grid      = element_blank()
    )
)

save_pub <- function(p, name, w_mm, h_mm) {
  w <- w_mm / 25.4; h <- h_mm / 25.4
  svglite::svglite(file.path(outdir, paste0(name, ".svg")), width = w, height = h)
  print(p); dev.off()
  grDevices::cairo_pdf(file.path(outdir, paste0(name, ".pdf")), width = w, height = h, family = fam)
  print(p); dev.off()
  ragg::agg_tiff(file.path(outdir, paste0(name, ".tiff")), width = w, height = h, units = "in", res = 600)
  print(p); dev.off()
  ragg::agg_png(file.path(outdir, paste0(name, ".png")), width = w, height = h, units = "in", res = 300)
  print(p); dev.off()
}

# Fig A: TLS status by response (n = 19)
df1 <- data.frame(
  group = factor(rep(c("Responder", "Non-responder"), times = c(15, 4)),
                 levels = c("Responder", "Non-responder")),
  tls   = factor(c(rep("TLS+", 14), "TLS-", rep("TLS+", 2), rep("TLS-", 2)),
                 levels = c("TLS+", "TLS-"))
)

## Barnard's exact test (two-sided, |Wald| ordering; Barnard 1947)
## inline implementation: CRAN Barnard/Exact packages are archived
barnard_test <- function(x1, n1, x2, n2, dp = 1e-4) {
  wald <- function(a, b) {
    p1 <- a / n1; p2 <- b / n2; pp <- (a + b) / (n1 + n2)
    if (pp == 0 || pp == 1) return(0)
    (p1 - p2) / sqrt(pp * (1 - pp) * (1 / n1 + 1 / n2))
  }
  w_obs <- wald(x1, x2)
  g <- expand.grid(a = 0:n1, b = 0:n2)
  g$w <- mapply(wald, g$a, g$b)
  ext <- abs(g$w) >= abs(w_obs) - 1e-12
  pis <- seq(dp, 1 - dp, by = dp)
  max(sapply(pis, function(pi)
    sum(dbinom(g$a[ext], n1, pi) * dbinom(g$b[ext], n2, pi))))
}
p_barnard <- barnard_test(14, 15, 2, 4)   # TLS+ in R / NR
p_lab <- sprintf("Barnard's exact p = %.3f", p_barnard)
cat(p_lab, "\n")

top1 <- data.frame(
  group = factor(c("Responder", "Non-responder"), levels = levels(df1$group)),
  y = c(15, 4),
  lab = c("93% TLS+", "50% TLS+")
)

p1 <- ggplot(df1, aes(group, fill = tls)) +
  geom_bar(width = 0.62, colour = "black", linewidth = 0.25) +
  geom_text(stat = "count", position = position_stack(vjust = 0.5),
            aes(label = after_stat(count), colour = tls),
            size = 2.2, family = fam) +
  geom_text(data = top1, aes(group, y, label = lab),
            vjust = -0.5, size = 2.1, family = fam, inherit.aes = FALSE) +
  annotate("text", x = 1.5, y = 18.6, label = p_lab, size = 2.1, family = fam) +
  scale_fill_manual(values = c("TLS+" = c_red, "TLS-" = c_grey)) +
  scale_colour_manual(values = c("TLS+" = "white", "TLS-" = "grey25"), guide = "none") +
  scale_x_discrete(labels = c("Responder\n(n = 15)", "Non-responder\n(n = 4)")) +
  scale_y_continuous(breaks = c(0, 5, 10, 15),
                     expand = expansion(mult = c(0, 0.16))) +
  labs(x = NULL, y = "Number of patients")

save_pub(p1, "figA_tls_by_response", w_mm = 60, h_mm = 62)

# Fig B: pre-treatment Hub+ vs Hub- -> TLS / response
dfB <- data.frame(
  endpoint = factor(rep(c("TLS after treatment", "Objective response"), each = 4),
                    levels = c("TLS after treatment", "Objective response")),
  hub    = factor(rep(rep(c("Hub+", "Hub-"), each = 2), 2), levels = c("Hub+", "Hub-")),
  status = factor(rep(c("Yes", "No"), 4), levels = c("Yes", "No")),
  n      = c(13, 1, 3, 2,   13, 1, 2, 3)
)
dfB$fillgrp <- factor(ifelse(dfB$status == "No", "No", as.character(dfB$endpoint)),
                      levels = c("TLS after treatment", "Objective response", "No"))

pctB <- data.frame(
  endpoint = factor(rep(c("TLS after treatment", "Objective response"), each = 2),
                    levels = levels(dfB$endpoint)),
  hub = factor(rep(c("Hub+", "Hub-"), 2), levels = c("Hub+", "Hub-")),
  lab = c("93%", "60%", "93%", "40%")
)
pvalB <- data.frame(
  endpoint = factor(c("TLS after treatment", "Objective response"),
                    levels = levels(dfB$endpoint)),
  lab = c("Barnard p = 0.093", "Barnard p = 0.024")
)
segB <- data.frame(
  endpoint = factor(rep(levels(dfB$endpoint), each = 3), levels = levels(dfB$endpoint)),
  x    = c(1, 1, 2, 1, 1, 2),
  xend = c(2, 1, 2, 2, 1, 2),
  y    = c(1.14, 1.11, 1.11, 1.14, 1.11, 1.11),
  yend = c(1.14, 1.14, 1.14, 1.14, 1.14, 1.14)
)

p2 <- ggplot(dfB, aes(hub, n, fill = fillgrp)) +
  geom_col(position = "fill", width = 0.62, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = n, colour = fillgrp),
            position = position_fill(vjust = 0.5), size = 2.2, family = fam) +
  geom_text(data = pctB, aes(hub, 1.04, label = lab),
            size = 2.1, family = fam, inherit.aes = FALSE) +
  geom_text(data = pvalB, aes(x = 1.5, y = 1.22, label = lab),
            size = 2.1, family = fam, inherit.aes = FALSE) +
  geom_segment(data = segB, aes(x = x, xend = xend, y = y, yend = yend),
               linewidth = 0.3, inherit.aes = FALSE) +
  facet_wrap(~endpoint, nrow = 1) +
  scale_fill_manual(values = c("TLS after treatment" = c_red,
                               "Objective response" = c_blue,
                               "No" = c_grey),
                    labels = c("TLS+", "Responder", "Not achieved")) +
  scale_colour_manual(values = c("TLS after treatment" = "white",
                                 "Objective response" = "white",
                                 "No" = "grey25"), guide = "none") +
  scale_x_discrete(labels = c("Hub+\n(n = 14)", "Hub-\n(n = 5)")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0%", "50%", "100%"),
                     expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1.3)) +
  labs(x = NULL, y = "Proportion of patients") +
  theme(strip.text = element_text(size = 6.2, face = "bold"),
        legend.position = "top")

save_pub(p2, "figB_hub_groups_outcomes", w_mm = 95, h_mm = 62)

# Fig C: Hub -> TLS -> response chain
library(ggalluvial)
dfS <- data.frame(
  hub  = factor(c("Hub+ (14)", "Hub+ (14)",
                  rep(c("Hub- (5)", "Hub- (5)"), each = 2)),
                levels = c("Hub+ (14)", "Hub- (5)")),
  tls  = factor(c("TLS+ (16)", "TLS- (3)",
                  "TLS+ (16)", "TLS+ (16)", "TLS- (3)", "TLS- (3)"),
                levels = c("TLS+ (16)", "TLS- (3)")),
  resp = factor(c("Responder (15)", "Non-responder (4)",
                  "Responder (15)", "Non-responder (4)",
                  "Responder (15)", "Non-responder (4)"),
                levels = c("Responder (15)", "Non-responder (4)")),
  hubfill = c("Hub+", "Hub+", "Hub-", "Hub-", "Hub-", "Hub-"),
  n = c(13, 1, 1, 2, 1, 1)
)

strata_tint <- c("Hub+ (14)" = "#DCE4F2", "Hub- (5)" = "#ECECEC",
                 "TLS+ (16)" = "#F2E5E4", "TLS- (3)" = "#ECECEC",
                 "Responder (15)" = "#DCE4F2", "Non-responder (4)" = "#ECECEC")

p3 <- ggplot(dfS, aes(y = n, axis1 = hub, axis2 = tls, axis3 = resp)) +
  geom_alluvium(aes(fill = hubfill), width = 1/4, alpha = 0.9,
                curve_type = "sigmoid") +
  scale_fill_manual(values = c("Hub+" = c_blue, "Hub-" = "#C4C4C4"), guide = "none") +
  ggnewscale::new_scale_fill() +
  geom_stratum(aes(fill = after_stat(stratum)), width = 1/4,
               colour = "grey30", linewidth = 0.3) +
  scale_fill_manual(values = strata_tint, guide = "none") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 2.0, family = fam) +
  scale_x_discrete(limits = c("hub", "tls", "resp"),
                   labels = c("Pre-treatment\nHub", "Post-treatment\nTLS", "Treatment\nresponse"),
                   expand = c(0.08, 0.08)) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 6.5, base_family = fam) +
  theme(axis.text.x = element_text(size = 6.5, face = "bold", colour = "black"),
        plot.margin = margin(2, 14, 2, 10, "mm"))

save_pub(p3, "figC_hub_tls_response_chain", w_mm = 105, h_mm = 70)

cat("done\n")
