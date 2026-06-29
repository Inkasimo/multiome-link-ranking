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
              help = "Full LinkPeaks baseline table; used only for diagnostics / compatibility"),
  make_option(c("--baseline-distance-file"), type = "character", default = file.path("results", "pbmc", "features", "pbmc_baseline_links_with_distance.csv"),
              help = "Full LinkPeaks baseline table with distance columns; used only for diagnostics / compatibility"),
  make_option(c("--score-mode"), type = "character", default = "full",
              help = "One of: linkpeaks, coactivity, coactivity_distance, coactivity_tf, distance_only, full, full_linkpeaks_anchored"),
  make_option(c("--lambda-distance"), type = "double", default = 0.30,
              help = "Strength of distance modifier [default %default]"),
  make_option(c("--alpha-tf"), type = "double", default = 0.50,
              help = "Strength of TF modifier [default %default]"),
  make_option(c("--ora-top-n"), type = "integer", default = 100,
              help = "Top N gene-ranked genes for ORA [default %default]"),
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

require_cols <- function(dt, cols, label) {
  missing <- setdiff(cols, names(dt))
  if (length(missing) > 0) {
    stop(label, " missing required columns: ", paste(missing, collapse = ", "))
  }
}

rescale01 <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) {
    out <- rep(0, length(x))
    out[is.na(x)] <- 0
    return(out)
  }
  out <- (x - rng[1]) / diff(rng)
  out[is.na(out)] <- 0
  out
}

rank_desc <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (any(ok)) out[ok] <- frank(-x[ok], ties.method = "average")
  out
}

assign_tiers <- function(score_vec, high_q = 0.90, medium_q = 0.70) {
  if (medium_q >= high_q) stop("tier-medium-quantile must be < tier-high-quantile")

  score_vec <- as.numeric(score_vec)
  finite <- is.finite(score_vec)
  out <- rep("Low", length(score_vec))

  if (sum(finite) < 3 || length(unique(score_vec[finite])) < 3) {
    return(factor(out, levels = c("High", "Medium", "Low")))
  }

  hi <- as.numeric(stats::quantile(score_vec[finite], probs = high_q, na.rm = TRUE))
  med <- as.numeric(stats::quantile(score_vec[finite], probs = medium_q, na.rm = TRUE))

  out[finite & score_vec >= med] <- "Medium"
  out[finite & score_vec >= hi] <- "High"
  factor(out, levels = c("High", "Medium", "Low"))
}

safe_bitr <- function(genes) {
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (length(genes) == 0) {
    return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
  }

  unique(suppressWarnings(
    bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  ))
}

safe_enrich_go <- function(entrez_ids, universe_ids) {
  entrez_ids <- unique(as.character(entrez_ids))
  universe_ids <- unique(as.character(universe_ids))
  entrez_ids <- entrez_ids[!is.na(entrez_ids) & nzchar(entrez_ids)]
  universe_ids <- universe_ids[!is.na(universe_ids) & nzchar(universe_ids)]

  if (length(entrez_ids) == 0 || length(universe_ids) == 0) return(NULL)

  tryCatch({
    enrichGO(
      gene = entrez_ids,
      universe = universe_ids,
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
  if (is.null(ora_obj) || nrow(as.data.frame(ora_obj)) == 0) {
    return(empty_plot(title_txt))
  }

  dotplot(ora_obj, showCategory = show_n) +
    ggtitle(title_txt) +
    theme_bw()
}

top_pairs_by_score <- function(dt, score_col, n) {
  require_cols(dt, c("peak", "gene", score_col), "top_pairs_by_score input")
  tmp <- copy(dt)
  tmp[, score_tmp__ := as.numeric(get(score_col))]
  tmp <- tmp[is.finite(score_tmp__)]
  setorderv(tmp, "score_tmp__", order = -1)

  unique(head(tmp[, .(peak, gene)], n))
}

top_genes_by_score <- function(dt, score_col, n) {
  require_cols(dt, c("gene", score_col), "top_genes_by_score input")
  gene_rank <- make_gene_rank(dt, score_col)
  head(gene_rank$gene, n)
}

make_gene_rank <- function(dt, score_col) {
  require_cols(dt, c("gene", score_col), "make_gene_rank input")

  tmp <- copy(dt)
  tmp[, score_tmp__ := as.numeric(get(score_col))]
  tmp <- tmp[is.finite(score_tmp__)]

  if (nrow(tmp) == 0) {
    return(data.table(
      gene = character(0),
      gene_score = numeric(0),
      n_links = integer(0),
      best_peak = character(0),
      gene_rank = integer(0)
    ))
  }

  out <- tmp[, {
    best_i <- which.max(score_tmp__)
    .(
      gene_score = score_tmp__[best_i],
      n_links = .N,
      best_peak = peak[best_i]
    )
  }, by = gene]

  setorder(out, -gene_score, gene)
  out[, gene_rank := .I]
  out[]
}

coerce_feature_aliases <- function(features) {
  features <- copy(features)

  if (!"tf_score" %in% names(features) && "peak_tf_score" %in% names(features)) {
    features[, tf_score := peak_tf_score]
  }

  if (!"motif_score" %in% names(features) && "peak_motif_score" %in% names(features)) {
    features[, motif_score := peak_motif_score]
  }

  if (!"motif_names" %in% names(features) && "peak_top_motifs" %in% names(features)) {
    features[, motif_names := peak_top_motifs]
  }

  features
}

get_rank_table <- function(features, mode, lambda_distance, alpha_tf, high_q, medium_q) {
  # SAME_CANDIDATE_UNIVERSE:
  # Every score mode ranks the same feature-table peak-gene pairs.
  dt <- copy(features)

  if (mode == "linkpeaks") {
    dt[, model_score := link_score]

  } else if (mode == "distance_only") {
    dt[, model_score := distance_score]

  } else if (mode == "coactivity") {
    dt[, model_score := mul_weigh]

  } else if (mode == "coactivity_distance") {
    dt[, model_score :=
      mul_weigh *
      ((1 - lambda_distance) + lambda_distance * distance_score)
    ]

  } else if (mode == "coactivity_tf") {
    dt[, model_score :=
      mul_weigh *
      (1 + alpha_tf * tf_score)
    ]

  } else if (mode == "full") {
    dt[, model_score :=
      mul_weigh *
      ((1 - lambda_distance) + lambda_distance * distance_score) *
      (1 + alpha_tf * tf_score)
    ]

  } else if (mode == "full_linkpeaks_anchored") {
    dt[, link_score_scaled := rescale01(link_score)]
    dt[, mul_weigh_scaled := rescale01(mul_weigh)]
    dt[, model_score :=
      link_score_scaled *
      (0.5 + 0.5 * mul_weigh_scaled) *
      ((1 - lambda_distance) + lambda_distance * distance_score) *
      (1 + alpha_tf * tf_score)
    ]

  } else {
    stop("Unknown score-mode: ", mode)
  }

  dt[, score_mode := mode]
  dt[, lambda_distance := lambda_distance]
  dt[, alpha_tf := alpha_tf]

  dt[, rank_model := rank_desc(model_score)]
  dt[, rank_link := rank_desc(link_score)]
  dt[, rank_diff_vs_linkpeaks := rank_link - rank_model]
  dt[, tier := assign_tiers(model_score, high_q, medium_q)]

  dt[, rank_order_tmp__ := ifelse(is.na(rank_model), Inf, rank_model)]
  setorder(dt, rank_order_tmp__, peak, gene)
  dt[, rank_order_tmp__ := NULL]

  dt[]
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

features <- coerce_feature_aliases(features)

require_cols(features, c(
  "peak", "gene", "link_score", "mul_weigh",
  "distance_score", "distance_bp", "tf_score"
), "features-file")

require_cols(baseline_full, c("peak", "gene", "score"), "baseline-file")
require_cols(baseline_dist, c("peak", "gene", "score", "distance_score", "distance_bp"), "baseline-distance-file")

features <- features[!is.na(peak) & !is.na(gene)]
features <- features[is.finite(link_score)]

setorder(features, peak, gene, -link_score)
features <- features[, .SD[1], by = .(peak, gene)]
setorder(features, -link_score)

if (nrow(features) == 0) {
  stop("No eligible feature rows after filtering.")
}

setorder(baseline_full, -score)
setorder(baseline_dist, -score)

msg("Evaluating score mode: %s", opt$score_mode)
msg("Evaluation universe: feature-table candidate pairs only")
msg("Eligible candidate pairs: %d", nrow(features))
msg("Eligible candidate genes: %d", uniqueN(features$gene))

ranked <- get_rank_table(
  features = features,
  mode = opt$score_mode,
  lambda_distance = opt$lambda_distance,
  alpha_tf = opt$alpha_tf,
  high_q = opt$tier_high_quantile,
  medium_q = opt$tier_medium_quantile
)

ranked_file <- file.path(opt$output_dir, sprintf("%s_ranked_links.csv", opt$run_name))
fwrite(ranked, ranked_file)
msg("Wrote ranked links: %s", ranked_file)

# Restrict full LinkPeaks files to the same candidate universe for diagnostics.
ranked_pairs <- unique(ranked[, .(peak, gene)])

baseline_eval <- merge(
  baseline_full,
  ranked_pairs,
  by = c("peak", "gene")
)
setorder(baseline_eval, -score)

baseline_dist_eval <- merge(
  baseline_dist,
  ranked_pairs,
  by = c("peak", "gene")
)
setorder(baseline_dist_eval, -score)

fwrite(baseline_eval, file.path(opt$output_dir, sprintf("%s_linkpeaks_baseline_same_universe.csv", opt$run_name)))
fwrite(baseline_dist_eval, file.path(opt$output_dir, sprintf("%s_linkpeaks_baseline_distance_same_universe.csv", opt$run_name)))

# ============================================================
# Standard top tables
# ============================================================
top100 <- head(ranked, 100)
fwrite(top100, file.path(opt$output_dir, sprintf("%s_top100_links.csv", opt$run_name)))

promoted <- ranked[order(-rank_diff_vs_linkpeaks)]
demoted <- ranked[order(rank_diff_vs_linkpeaks)]

fwrite(head(promoted, 100), file.path(opt$output_dir, sprintf("%s_top_promoted_vs_linkpeaks.csv", opt$run_name)))
fwrite(head(demoted, 100), file.path(opt$output_dir, sprintf("%s_top_demoted_vs_linkpeaks.csv", opt$run_name)))

tier_summary <- ranked[, .(
  n_links = .N,
  unique_genes = uniqueN(gene),
  median_score = median(model_score, na.rm = TRUE),
  median_link_score = median(link_score, na.rm = TRUE),
  median_distance_bp = median(distance_bp[is.finite(distance_bp)], na.rm = TRUE)
), by = tier]
fwrite(tier_summary, file.path(opt$output_dir, sprintf("%s_tier_summary.csv", opt$run_name)))

model_gene_rank <- make_gene_rank(ranked, "model_score")
linkpeaks_gene_rank <- make_gene_rank(ranked, "link_score")

fwrite(model_gene_rank, file.path(opt$output_dir, sprintf("%s_model_gene_rank.csv", opt$run_name)))
fwrite(linkpeaks_gene_rank, file.path(opt$output_dir, sprintf("%s_linkpeaks_gene_rank_same_universe.csv", opt$run_name)))

# ============================================================
# ORA
# ============================================================
msg("Running ORA on top %d gene-ranked genes...", opt$ora_top_n)

model_genes <- head(model_gene_rank$gene, opt$ora_top_n)
baseline_genes <- head(linkpeaks_gene_rank$gene, opt$ora_top_n)
background_genes <- unique(ranked$gene)

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
fwrite(baseline_ora_df, file.path(opt$output_dir, sprintf("%s_linkpeaks_baseline_ora_GO_BP_same_universe.csv", opt$run_name)))

model_dot <- make_dotplot_safe(ora_model, sprintf("ORA: %s", opt$score_mode), opt$ora_show_category)
ggsave(file.path(opt$output_dir, sprintf("%s_ora_dotplot.png", opt$run_name)), model_dot, width = 10, height = 7, dpi = 300)

baseline_dot <- make_dotplot_safe(ora_baseline, "ORA: LinkPeaks baseline, same universe", opt$ora_show_category)
ggsave(file.path(opt$output_dir, sprintf("%s_linkpeaks_baseline_ora_dotplot_same_universe.png", opt$run_name)), baseline_dot, width = 10, height = 7, dpi = 300)

# ============================================================
# Metrics
# ============================================================
ranked_model_finite_dist <- ranked[is.finite(distance_bp)]
ranked_link_finite_dist <- copy(ranked_model_finite_dist)
setorder(ranked_link_finite_dist, -link_score)

top50_model <- head(ranked_model_finite_dist, 50)
top50_link <- head(ranked_link_finite_dist, 50)

top100_model <- top_pairs_by_score(ranked, "model_score", 100)
top100_link <- top_pairs_by_score(ranked, "link_score", 100)

pair_overlap <- nrow(merge(as.data.frame(top100_model), as.data.frame(top100_link), by = c("peak", "gene")))

top100_model_genes <- top_genes_by_score(ranked, "model_score", 100)
top100_link_genes <- top_genes_by_score(ranked, "link_score", 100)
gene_overlap <- length(intersect(top100_model_genes, top100_link_genes))

summary_dt <- data.table(
  metric = c(
    "score_mode",
    "evaluation_universe",
    "n_links_ranked",
    "n_unique_genes_ranked",
    "n_full_linkpeaks_links_unrestricted",
    "n_linkpeaks_links_same_universe",
    "cor_linkpeaks_vs_model_score_same_universe",
    "median_distance_top50_linkpeaks_same_universe",
    "median_distance_top50_model",
    "distal_frac_top50_linkpeaks_gt50kb_same_universe",
    "distal_frac_top50_model_gt50kb",
    "n_unique_linkpeaks_genes_topN_gene_ranked",
    "n_unique_model_genes_topN_gene_ranked",
    "n_ora_background_genes",
    "n_baseline_ora_terms_same_universe",
    "n_model_ora_terms",
    "overlap_pairs_top100_model_vs_linkpeaks_same_universe",
    "overlap_genes_top100_model_vs_linkpeaks_same_universe"
  ),
  value = c(
    opt$score_mode,
    "features_same_candidate_universe",
    as.character(nrow(ranked)),
    as.character(uniqueN(ranked$gene)),
    as.character(nrow(baseline_full)),
    as.character(nrow(baseline_eval)),
    as.character(suppressWarnings(cor(ranked$link_score, ranked$model_score, use = "complete.obs"))),
    as.character(median(top50_link$distance_bp, na.rm = TRUE)),
    as.character(median(top50_model$distance_bp, na.rm = TRUE)),
    as.character(mean(top50_link$distance_bp > 50000, na.rm = TRUE)),
    as.character(mean(top50_model$distance_bp > 50000, na.rm = TRUE)),
    as.character(length(unique(baseline_genes))),
    as.character(length(unique(model_genes))),
    as.character(length(unique(background_genes))),
    as.character(nrow(baseline_ora_df)),
    as.character(nrow(model_ora_df)),
    as.character(pair_overlap),
    as.character(gene_overlap)
  )
)

fwrite(summary_dt, file.path(opt$output_dir, sprintf("%s_summary_metrics.csv", opt$run_name)))

# Top-N overlap summary within the same candidate universe.
noverlap <- rbindlist(lapply(c(10, 20, 50, 100, 200), function(n) {
  top_link_n <- top_pairs_by_score(ranked, "link_score", n)
  top_model_n <- top_pairs_by_score(ranked, "model_score", n)
  overlap_n <- nrow(merge(as.data.frame(top_link_n), as.data.frame(top_model_n), by = c("peak", "gene")))

  data.table(top_n = n, overlap = overlap_n)
}))

fwrite(noverlap, file.path(opt$output_dir, sprintf("%s_topN_overlap_vs_linkpeaks_same_universe.csv", opt$run_name)))

# ============================================================
# Plots
# ============================================================
p_scatter <- ggplot(ranked, aes(x = link_score, y = model_score, color = tier)) +
  geom_point(alpha = 0.5, size = 1) +
  theme_bw() +
  ggtitle(sprintf("LinkPeaks vs %s, same candidate universe", opt$score_mode)) +
  xlab("LinkPeaks score") +
  ylab(sprintf("%s score", opt$score_mode))

ggsave(file.path(opt$output_dir, sprintf("%s_linkpeaks_vs_model_scatter.png", opt$run_name)), p_scatter, width = 8, height = 6, dpi = 300)

ranked[, top200_model := rank_model <= 200]

p_dist <- ggplot(ranked[is.finite(distance_bp)], aes(x = distance_bp, fill = top200_model)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  theme_bw() +
  ggtitle(sprintf("Distance distribution: all links vs top200 %s", opt$score_mode)) +
  xlab("Distance to nearest transcript TSS (bp)") +
  ylab("Count")

ggsave(file.path(opt$output_dir, sprintf("%s_distance_distribution.png", opt$run_name)), p_dist, width = 8, height = 6, dpi = 300)

if (!is.null(opt$done_file)) {
  dir.create(dirname(opt$done_file), showWarnings = FALSE, recursive = TRUE)
  file.create(opt$done_file)
  msg("Touched done file: %s", opt$done_file)
}

msg("Done.")
msg("Main output: %s", ranked_file)
