#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
})

# ============================================================
# CLI
# ============================================================
option_list <- list(
  make_option(c("--features-file"), type = "character", default = file.path("results", "pbmc", "features", "pbmc_link_features.csv"),
              help = "Feature table produced by run_linkpeaks_reranker.R"),
  make_option(c("--baseline-file"), type = "character", default = file.path("results", "pbmc", "features", "pbmc_baseline_links_full.csv"),
              help = "Full LinkPeaks baseline table"),
  make_option(c("--baseline-distance-file"), type = "character", default = file.path("results", "pbmc", "features", "pbmc_baseline_links_with_distance.csv"),
              help = "Full LinkPeaks baseline table with distance columns"),
  make_option(c("--score-mode"), type = "character", default = "full",
              help = "One of: linkpeaks, coactivity, coactivity_distance, coactivity_tf, distance_only, full"),
  make_option(c("--lambda-distance"), type = "double", default = 0.30,
              help = "Strength of distance modifier [default %default]"),
  make_option(c("--alpha-tf"), type = "double", default = 0.50,
              help = "Strength of TF modifier [default %default]"),
  make_option(c("--ora-top-n"), type = "integer", default = 100,
              help = "Top N genes for ORA [default %default]"),
  make_option(c("--ora-show-category"), type = "integer", default = 15,
              help = "Number of categories to show in ORA dotplots [default %default]"),
  make_option(c("--tier-high-quantile"), type = "double", default = 0.90,
              help = "Quantile cutoff for High tier [default %default]"),
  make_option(c("--tier-medium-quantile"), type = "double", default = 0.70,
              help = "Quantile cutoff for Medium tier [default %default]"),
  make_option(c("--output-dir"), type = "character", default = file.path("results", "pbmc", "rankings", "full"),
              help = "Directory for ranking/evaluation outputs"),
  make_option(c("--run-name"), type = "character", default = "pbmc_full",
              help = "Prefix for output files [default %default]"),
  make_option(c("--done-file"), type = "character", default = NULL,
              help = "Optional sentinel file touched on success")
)

opt <- parse_args(OptionParser(option_list = option_list))
names(opt) <- gsub("-", "_", names(opt), fixed = TRUE)

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Helpers
# ============================================================
msg <- function(...) cat(sprintf(...), "\n")

require_file <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
}

assign_tiers <- function(score_vec, high_q = 0.90, medium_q = 0.70) {
  if (medium_q >= high_q) stop("tier-medium-quantile must be < tier-high-quantile")
  hi <- as.numeric(stats::quantile(score_vec, probs = high_q, na.rm = TRUE))
  med <- as.numeric(stats::quantile(score_vec, probs = medium_q, na.rm = TRUE))
  out <- ifelse(score_vec >= hi, "High", ifelse(score_vec >= med, "Medium", "Low"))
  factor(out, levels = c("High", "Medium", "Low"))
}

safe_bitr <- function(genes) {
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (length(genes) == 0) return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
  unique(suppressWarnings(
    bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  ))
}

safe_enrich_go <- function(entrez_ids, universe_ids) {
  if (length(entrez_ids) == 0 || length(universe_ids) == 0) return(NULL)
  tryCatch({
    enrichGO(
      gene = unique(entrez_ids),
      universe = unique(universe_ids),
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.20,
      readable = TRUE
    )
  }, error = function(e) {
    warning("enrichGO failed: ", conditionMessage(e))
    NULL
  })
}

empty_plot <- function(title_txt) {
  ggplot() +
    annotate("text", x = 1, y = 1, label = "No significant enrichment", size = 5) +
    xlim(0, 2) + ylim(0, 2) +
    theme_void() +
    ggtitle(title_txt)
}

make_dotplot_safe <- function(ora_obj, title_txt, show_n = 15) {
  if (is.null(ora_obj) || nrow(as.data.frame(ora_obj)) == 0) return(empty_plot(title_txt))
  dotplot(ora_obj, showCategory = show_n) + ggtitle(title_txt) + theme_bw()
}

pair_id <- function(dt) paste(dt$peak, dt$gene, sep = "||")

get_rank_table <- function(features, baseline_dist, mode, lambda_distance, alpha_tf) {
  if (mode == "linkpeaks") {
    dt <- copy(baseline_dist)
    setnames(dt, "score", "link_score", skip_absent = TRUE)
    dt[, model_score := link_score]
  } else if (mode == "distance_only") {
    dt <- copy(baseline_dist)
    setnames(dt, "score", "link_score", skip_absent = TRUE)
    dt[, model_score := distance_score]
  } else {
    dt <- copy(features)
    if (mode == "coactivity") {
      dt[, model_score := mul_weigh]
    } else if (mode == "coactivity_distance") {
      dt[, model_score := mul_weigh * ((1 - lambda_distance) + lambda_distance * distance_score)]
    } else if (mode == "coactivity_tf") {
      dt[, model_score := mul_weigh * (1 + alpha_tf * tf_score)]
    } else if (mode == "full") {
      dt[, model_score := mul_weigh * ((1 - lambda_distance) + lambda_distance * distance_score) * (1 + alpha_tf * tf_score)]
    } else {
      stop("Unknown score-mode: ", mode)
    }
  }

  dt[, score_mode := mode]
  dt[, lambda_distance := lambda_distance]
  dt[, alpha_tf := alpha_tf]
  dt[, rank_model := frank(-model_score, ties.method = "average")]
  if ("link_score" %in% names(dt)) {
    dt[, rank_link := frank(-link_score, ties.method = "average")]
    dt[, rank_diff_vs_linkpeaks := rank_link - rank_model]
  }
  dt[, tier := assign_tiers(model_score, opt$tier_high_quantile, opt$tier_medium_quantile)]
  setorder(dt, rank_model)
  dt
}

# ============================================================
# Load inputs
# ============================================================
require_file(opt$features_file)
require_file(opt$baseline_file)
require_file(opt$baseline_distance_file)

features <- fread(opt$features_file)
baseline_full <- fread(opt$baseline_file)
baseline_dist <- fread(opt$baseline_distance_file)

if (!"distance_score" %in% names(baseline_dist)) {
  stop("baseline-distance-file must contain distance_score. Re-run run_linkpeaks_reranker.R.")
}

msg("Evaluating score mode: %s", opt$score_mode)
ranked <- get_rank_table(features, baseline_dist, opt$score_mode, opt$lambda_distance, opt$alpha_tf)

ranked_file <- file.path(opt$output_dir, sprintf("%s_ranked_links.csv", opt$run_name))
fwrite(ranked, ranked_file)
msg("Wrote ranked links: %s", ranked_file)

# ============================================================
# Standard top tables
# ============================================================
top100 <- head(ranked, 100)
fwrite(top100, file.path(opt$output_dir, sprintf("%s_top100_links.csv", opt$run_name)))

if ("rank_diff_vs_linkpeaks" %in% names(ranked)) {
  promoted <- ranked[order(-rank_diff_vs_linkpeaks)]
  demoted <- ranked[order(rank_diff_vs_linkpeaks)]
  fwrite(head(promoted, 100), file.path(opt$output_dir, sprintf("%s_top_promoted_vs_linkpeaks.csv", opt$run_name)))
  fwrite(head(demoted, 100), file.path(opt$output_dir, sprintf("%s_top_demoted_vs_linkpeaks.csv", opt$run_name)))
}

tier_summary <- ranked[, .(
  n_links = .N,
  unique_genes = uniqueN(gene),
  median_score = median(model_score, na.rm = TRUE),
  median_distance_bp = median(distance_bp[is.finite(distance_bp)], na.rm = TRUE)
), by = tier]
fwrite(tier_summary, file.path(opt$output_dir, sprintf("%s_tier_summary.csv", opt$run_name)))

# ============================================================
# ORA
# ============================================================
msg("Running ORA on top %d genes...", opt$ora_top_n)
model_genes <- unique(head(ranked$gene, opt$ora_top_n))
baseline_genes <- unique(head(baseline_full$gene, opt$ora_top_n))
background_genes <- unique(baseline_full$gene)

model_map <- safe_bitr(model_genes)
baseline_map <- safe_bitr(baseline_genes)
background_map <- safe_bitr(background_genes)

model_entrez <- unique(model_map$ENTREZID)
baseline_entrez <- unique(baseline_map$ENTREZID)
background_entrez <- unique(background_map$ENTREZID)

ora_model <- safe_enrich_go(model_entrez, background_entrez)
ora_baseline <- safe_enrich_go(baseline_entrez, background_entrez)

model_ora_df <- if (is.null(ora_model)) data.frame() else as.data.frame(ora_model)
baseline_ora_df <- if (is.null(ora_baseline)) data.frame() else as.data.frame(ora_baseline)

fwrite(model_ora_df, file.path(opt$output_dir, sprintf("%s_ora_GO_BP.csv", opt$run_name)))
fwrite(baseline_ora_df, file.path(opt$output_dir, sprintf("%s_baseline_ora_GO_BP.csv", opt$run_name)))

model_dot <- make_dotplot_safe(ora_model, sprintf("ORA: %s", opt$score_mode), opt$ora_show_category)
ggsave(file.path(opt$output_dir, sprintf("%s_ora_dotplot.png", opt$run_name)), model_dot, width = 10, height = 7, dpi = 300)

# ============================================================
# Metrics
# ============================================================
baseline_dist_finite <- baseline_dist[is.finite(distance_bp)]
setorder(baseline_dist_finite, -score)
baseline_top50 <- head(baseline_dist_finite, 50)

ranked_finite <- ranked[is.finite(distance_bp)]
top50 <- head(ranked_finite, 50)
top100_model <- unique(head(ranked[, .(peak, gene)], 100))
top100_link <- unique(head(baseline_full[, .(peak, gene)], 100))

pair_overlap <- nrow(merge(as.data.frame(top100_model), as.data.frame(top100_link), by = c("peak", "gene")))
gene_overlap <- length(intersect(unique(head(ranked$gene, 100)), unique(head(baseline_full$gene, 100))))

summary_dt <- data.table(
  metric = c(
    "score_mode",
    "n_links_ranked",
    "cor_linkpeaks_vs_model_score",
    "median_distance_top50_linkpeaks",
    "median_distance_top50_model",
    "distal_frac_top50_model_gt50kb",
    "n_unique_linkpeaks_genes_topN",
    "n_unique_model_genes_topN",
    "n_baseline_ora_terms",
    "n_model_ora_terms",
    "overlap_pairs_top100_model_vs_linkpeaks",
    "overlap_genes_top100_model_vs_linkpeaks"
  ),
  value = c(
    opt$score_mode,
    as.character(nrow(ranked)),
    as.character(suppressWarnings(cor(ranked$link_score, ranked$model_score, use = "complete.obs"))),
    as.character(median(baseline_top50$distance_bp, na.rm = TRUE)),
    as.character(median(top50$distance_bp, na.rm = TRUE)),
    as.character(mean(top50$distance_bp > 50000, na.rm = TRUE)),
    as.character(length(unique(baseline_genes))),
    as.character(length(unique(model_genes))),
    as.character(nrow(baseline_ora_df)),
    as.character(nrow(model_ora_df)),
    as.character(pair_overlap),
    as.character(gene_overlap)
  )
)
fwrite(summary_dt, file.path(opt$output_dir, sprintf("%s_summary_metrics.csv", opt$run_name)))

# Top-N overlap summary
noverlap <- rbindlist(lapply(c(10, 20, 50, 100, 200), function(n) {
  top_link_n <- unique(head(baseline_full[, .(peak, gene)], n))
  top_model_n <- unique(head(ranked[, .(peak, gene)], n))
  overlap_n <- nrow(merge(as.data.frame(top_link_n), as.data.frame(top_model_n), by = c("peak", "gene")))
  data.table(top_n = n, overlap = overlap_n)
}))
fwrite(noverlap, file.path(opt$output_dir, sprintf("%s_topN_overlap_vs_linkpeaks.csv", opt$run_name)))

# ============================================================
# Plots
# ============================================================
if ("link_score" %in% names(ranked)) {
  p_scatter <- ggplot(ranked, aes(x = link_score, y = model_score, color = tier)) +
    geom_point(alpha = 0.5, size = 1) +
    theme_bw() +
    ggtitle(sprintf("LinkPeaks vs %s", opt$score_mode)) +
    xlab("LinkPeaks score") +
    ylab(sprintf("%s score", opt$score_mode))
  ggsave(file.path(opt$output_dir, sprintf("%s_linkpeaks_vs_model_scatter.png", opt$run_name)), p_scatter, width = 8, height = 6, dpi = 300)
}

ranked[, top200_model := seq_len(.N) <= 200]
p_dist <- ggplot(ranked[is.finite(distance_bp)], aes(x = distance_bp, fill = top200_model)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  theme_bw() +
  ggtitle(sprintf("Distance distribution: all links vs top200 %s", opt$score_mode)) +
  xlab("Distance to gene TSS (bp)") +
  ylab("Count")
ggsave(file.path(opt$output_dir, sprintf("%s_distance_distribution.png", opt$run_name)), p_dist, width = 8, height = 6, dpi = 300)

if (!is.null(opt$done_file)) {
  dir.create(dirname(opt$done_file), showWarnings = FALSE, recursive = TRUE)
  file.create(opt$done_file)
  msg("Touched done file: %s", opt$done_file)
}

msg("Done.")
msg("Main output: %s", ranked_file)
