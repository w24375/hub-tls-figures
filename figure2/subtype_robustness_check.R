#!/usr/bin/env Rscript
# Robustness check for Fig 2: reclassify Base/Stem hubs with TLS-signature
# genes REMOVED from the Hub_stem module, then recompute the sample-level
# distance comparison.
#
# Original Hub_stem = TCF7, PDCD1, CCR7, CCL19, CCL17, CCL22.
# TLS detection modules contribute CCL19 (Sautes_9/Chemo_core/E-TLS), CCR7
# (E-TLS/detection), PDCD1 (Tcell_mod/PFL-TLS).
#   strict3 : TCF7 + CCL17 + CCL22          (no overlap with any TLS module)
#   nocc4   : TCF7 + PDCD1 + CCL17 + CCL22  (only E-TLS overlap removed)
# Classification semantics identical to run_full_tls_pipeline.R:
#   stem = is_immune_hub & score >= median(score, ALL spots in sample).

suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2)
})

set.seed(20261007)
OUT <- "data/hub_tls_paper_outline/figure2_assembled/robustness_checks"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

dist_dt <- fread("data/spatial/hub_precursor_analysis/spatial_gradient/spot_spatial_distances.csv")
samples <- unique(dist_dt[, .(dataset, sample)])
cat(sprintf("Samples in distance table: %d\n", nrow(samples)))

score_one <- function(dataset, sample) {
  rds <- file.path("data/spatial/organized_rds", dataset, paste0(sample, ".Rds"))
  if (!file.exists(rds)) return(NULL)
  seu <- readRDS(rds)
  if (!"data" %in% slotNames(seu@assays[[DefaultAssay(seu)]])) {
    seu <- NormalizeData(seu, verbose = FALSE)
  }
  seu <- AddModuleScore(seu, features = list(c("TCF7", "CCL17", "CCL22")),
                        name = "strict3", ctrl = 15)
  seu <- AddModuleScore(seu, features = list(c("TCF7", "PDCD1", "CCL17", "CCL22")),
                        name = "nocc4", ctrl = 20)
  md <- seu@meta.data
  out <- data.table(
    barcode  = paste(sample, colnames(seu), sep = "_"),
    strict3  = md$strict31,
    nocc4    = md$nocc41,
    med_strict3 = median(md$strict31, na.rm = TRUE),
    med_nocc4   = median(md$nocc41,   na.rm = TRUE)
  )
  out
}

res <- list(); fail <- character()
for (i in seq_len(nrow(samples))) {
  d <- samples$dataset[i]; s <- samples$sample[i]
  r <- tryCatch(score_one(d, s), error = function(e) NULL)
  if (is.null(r)) { fail <- c(fail, paste(d, s, sep = "/")); next }
  res[[length(res) + 1]] <- r
  if (length(res) %% 25 == 0) cat(sprintf("  scored %d samples\n", length(res)))
}
scores <- rbindlist(res)
cat(sprintf("Scored %d samples (%d spots); failed: %d\n",
            length(res), nrow(scores), length(fail)))
if (length(fail)) cat("  failed:", paste(fail, collapse = ", "), "\n")

m <- merge(dist_dt, scores, by = "barcode", all.x = TRUE)
m <- m[!is.na(strict3)]
cat(sprintf("Merged: %d spots with strict scores\n", nrow(m)))

# Reclassify within original immune hubs (hub membership itself unchanged)
m[, stem_strict3 := is_immune_hub == TRUE & strict3 >= med_strict3]
m[, stem_nocc4   := is_immune_hub == TRUE & nocc4   >= med_nocc4]

recompute <- function(dt, stem_col) {
  h <- dt[is_immune_hub == TRUE]
  h[, sub := ifelse(get(stem_col), "Stem", "Base"), ]
  cnt <- h[, .(n_s = sum(sub == "Stem"), n_b = sum(sub == "Base")), by = .(sample, dataset)]
  keep <- cnt[n_s >= 5 & n_b >= 5, sample]
  h <- h[sample %in% keep]
  agg <- h[, .(mean_s = mean(dist_to_tls[sub == "Stem"]),
               mean_b = mean(dist_to_tls[sub == "Base"])), by = .(sample, dataset)]
  agg <- agg[mean_s > 0 & mean_b > 0]
  agg[, log2fc := log2(mean_s / mean_b)]
  agg
}

summarise_variant <- function(agg, tag) {
  data.table(
    variant = tag,
    n_samples = nrow(agg),
    pct_stem_closer = round(100 * mean(agg$log2fc < 0), 1),
    mean_log2fc = round(mean(agg$log2fc), 3),
    wilcox_p = signif(wilcox.test(agg$mean_s, agg$mean_b, paired = TRUE,
                                  alternative = "less")$p.value, 3),
    median_dist_stem = round(median(agg$mean_s), 1),
    median_dist_base = round(median(agg$mean_b), 1)
  )
}

agg_s3 <- recompute(m, "stem_strict3")
agg_n4 <- recompute(m, "stem_nocc4")

# Original labels for reference (same sample filter)
orig <- dist_dt[is_immune_hub == TRUE & !is.na(is_stem_hub)]
orig[, sub := ifelse(is_stem_hub == TRUE, "Stem", "Base")]
cnto <- orig[, .(n_s = sum(sub == "Stem"), n_b = sum(sub == "Base")), by = sample]
orig <- orig[sample %in% cnto[n_s >= 5 & n_b >= 5, sample]]
agg_o <- orig[, .(mean_s = mean(dist_to_tls[sub == "Stem"]),
                  mean_b = mean(dist_to_tls[sub == "Base"])), by = sample]
agg_o <- agg_o[mean_s > 0 & mean_b > 0][, log2fc := log2(mean_s / mean_b)]

stats <- rbind(summarise_variant(agg_o,  "original (Hub_stem 6 genes)"),
               summarise_variant(agg_n4, "nocc4 (CCL19/CCR7 removed)"),
               summarise_variant(agg_s3, "strict3 (CCL19/CCR7/PDCD1 removed)"))
print(stats)
fwrite(stats, file.path(OUT, "subtype_robustness_stats.csv"))
fwrite(agg_o,  file.path(OUT, "per_sample_original.csv"))
fwrite(agg_n4, file.path(OUT, "per_sample_nocc4.csv"))
fwrite(agg_s3, file.path(OUT, "per_sample_strict3.csv"))

# Figure: per-sample log2FC under original vs de-TLS classification
mm <- merge(agg_o[, .(sample, orig = log2fc)], agg_s3[, .(sample, strict3 = log2fc)], by = "sample")
gp <- ggplot(mm, aes(x = orig, y = strict3)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_point(size = 1.2, alpha = 0.7, colour = "#6C8EBF") +
  labs(x = "log2FC Stem/Base distance (original classification)",
       y = "log2FC Stem/Base distance\n(TLS-gene-free classification)",
       title = "Stem-closer-to-TLS direction survives removal of TLS-signature genes") +
  theme_classic(base_size = 8) +
  theme(plot.title = element_text(size = 8, face = "bold"))
ggsave(file.path(OUT, "subtype_robustness_scatter.pdf"), gp, width = 95, height = 90, units = "mm")
ggsave(file.path(OUT, "subtype_robustness_scatter.png"), gp, width = 95, height = 90, units = "mm", dpi = 300)
cat("Done.\n")
