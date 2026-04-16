#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# =========================
# Config
# =========================
run_dirs <- c("50000_03", "75000_03", "100000_03")

ranked_file  <- "multiome_rie_ranked_links.csv"
summary_file <- "multiome_rie_summary_metrics.csv"

top_n_main <- 100
top_n_alt  <- 200
high_dist_bp <- 50000
mid_dist_min <- 10000
mid_dist_max <- 50000
high_rank_cutoff <- 500

# =========================
# Helpers
# =========================
require_file <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
}

safe_read <- function(path) {
  require_file(path)
  fread(path)
}

metric_value <- function(summary_dt, metric_name) {
  x <- summary_dt[metric == metric_name, value]
  if (length(x) == 0) return(NA_real_)
  as.numeric(x[1])
}

pair_id <- function(dt) {
  paste(dt$peak, dt$gene, sep = "||")
}

top_n_dt <- function(dt, n = 100) {
  dt[order(rank_final_v6)][1:min(.N, n)]
}

high_dist_high_rank_dt <- function(dt) {
  dt[distance_bp > high_dist_bp & rank_final_v6 <= high_rank_cutoff][order(rank_final_v6)]
}

mid_dist_high_rank_dt <- function(dt) {
  dt[distance_bp > mid_dist_min & distance_bp < mid_dist_max & rank_final_v6 <= high_rank_cutoff][order(rank_final_v6)]
}

summarize_run <- function(run_name, ranked_dt, summary_dt) {
  top50  <- top_n_dt(ranked_dt, 50)
  top100 <- top_n_dt(ranked_dt, 100)
  top200 <- top_n_dt(ranked_dt, 200)

  high_dist <- high_dist_high_rank_dt(ranked_dt)
  mid_dist  <- mid_dist_high_rank_dt(ranked_dt)

  data.table(
    run = run_name,
    n_links = nrow(ranked_dt),

    median_top50_dist = median(top50$distance_bp, na.rm = TRUE),
    median_top100_dist = median(top100$distance_bp, na.rm = TRUE),

    distal_frac_top50_gt50kb = mean(top50$distance_bp > high_dist_bp, na.rm = TRUE),
    distal_frac_top100_gt50kb = mean(top100$distance_bp > high_dist_bp, na.rm = TRUE),
    promoter_frac_top100_lt2kb = mean(top100$distance_bp < 2000, na.rm = TRUE),

    unique_genes_top100 = uniqueN(top100$gene),
    unique_genes_top200 = uniqueN(top200$gene),

    median_tf_top100 = median(top100$tf_score, na.rm = TRUE),
    median_mul_weigh_top100 = median(top100$mul_weigh, na.rm = TRUE),

    n_high_dist_high_rank = nrow(high_dist),
    median_dist_high_dist_high_rank = if (nrow(high_dist) > 0) median(high_dist$distance_bp, na.rm = TRUE) else NA_real_,
    median_tf_high_dist_high_rank = if (nrow(high_dist) > 0) median(high_dist$tf_score, na.rm = TRUE) else NA_real_,
    median_mul_weigh_high_dist_high_rank = if (nrow(high_dist) > 0) median(high_dist$mul_weigh, na.rm = TRUE) else NA_real_,

    n_mid_dist_high_rank = nrow(mid_dist),
    median_tf_mid_dist_high_rank = if (nrow(mid_dist) > 0) median(mid_dist$tf_score, na.rm = TRUE) else NA_real_,
    median_mul_weigh_mid_dist_high_rank = if (nrow(mid_dist) > 0) median(mid_dist$mul_weigh, na.rm = TRUE) else NA_real_,

    cor_link_vs_final_v6 = metric_value(summary_dt, "cor_link_vs_final_v6"),
    n_unique_baseline_genes_topN = metric_value(summary_dt, "n_unique_baseline_genes_topN"),
    n_unique_final_genes_topN = metric_value(summary_dt, "n_unique_final_genes_topN"),
    n_baseline_ora_terms = metric_value(summary_dt, "n_baseline_ora_terms"),
    n_final_ora_terms = metric_value(summary_dt, "n_final_ora_terms")
  )
}

compare_top_sets <- function(run_data, n = 100) {
  out <- list()
  idx <- 1

  top_pairs <- list()
  top_genes <- list()

  for (nm in names(run_data)) {
    dt <- top_n_dt(run_data[[nm]]$ranked, n)
    top_pairs[[nm]] <- unique(pair_id(dt))
    top_genes[[nm]] <- unique(dt$gene)
  }

  runs <- names(run_data)
  for (i in seq_along(runs)) {
    for (j in seq_along(runs)) {
      if (j <= i) next
      a <- runs[i]
      b <- runs[j]

      pair_overlap <- length(intersect(top_pairs[[a]], top_pairs[[b]]))
      pair_union   <- length(union(top_pairs[[a]], top_pairs[[b]]))

      gene_overlap <- length(intersect(top_genes[[a]], top_genes[[b]]))
      gene_union   <- length(union(top_genes[[a]], top_genes[[b]]))

      out[[idx]] <- data.table(
        run_a = a,
        run_b = b,
        top_n = n,
        pair_overlap = pair_overlap,
        pair_jaccard = ifelse(pair_union > 0, pair_overlap / pair_union, NA_real_),
        gene_overlap = gene_overlap,
        gene_jaccard = ifelse(gene_union > 0, gene_overlap / gene_union, NA_real_)
      )
      idx <- idx + 1
    }
  }

  rbindlist(out)
}

top_examples <- function(dt, n = 20) {
  keep <- c("gene", "peak", "link_score", "mul_weigh", "distance_bp",
            "distance_score", "tf_score", "final_v6", "rank_link",
            "rank_final_v6", "rank_diff_v6", "tier", "motif_names")
  top_n_dt(dt, n)[, ..keep]
}

# =========================
# Load runs
# =========================
run_data <- list()
summary_rows <- list()

for (run_name in run_dirs) {
  ranked_path  <- file.path(run_name, ranked_file)
  summary_path <- file.path(run_name, summary_file)

  ranked_dt  <- safe_read(ranked_path)
  summary_dt <- safe_read(summary_path)

  run_data[[run_name]] <- list(
    ranked = ranked_dt,
    summary = summary_dt
  )

  summary_rows[[run_name]] <- summarize_run(run_name, ranked_dt, summary_dt)

  fwrite(
    top_examples(high_dist_high_rank_dt(ranked_dt), 20),
    file = sprintf("comparison_%s_high_distance_high_rank_links.csv", run_name)
  )

  fwrite(
    top_examples(mid_dist_high_rank_dt(ranked_dt), 20),
    file = sprintf("comparison_%s_mid_distance_high_rank_links.csv", run_name)
  )

  fwrite(
    top_examples(ranked_dt, 20),
    file = sprintf("comparison_%s_top20_final_links.csv", run_name)
  )
}

# =========================
# Main outputs
# =========================
comparison_summary <- rbindlist(summary_rows)
setorder(comparison_summary, run)

top100_overlap <- compare_top_sets(run_data, n = top_n_main)
top200_overlap <- compare_top_sets(run_data, n = top_n_alt)

fwrite(comparison_summary, "comparison_run_summary.csv")
fwrite(top100_overlap, "comparison_top100_overlap_between_runs.csv")
fwrite(top200_overlap, "comparison_top200_overlap_between_runs.csv")

# =========================
# Console output
# =========================
cat("\n=== RUN SUMMARY ===\n")
print(comparison_summary)

cat("\n=== TOP100 OVERLAP BETWEEN RUNS ===\n")
print(top100_overlap)

cat("\n=== TOP200 OVERLAP BETWEEN RUNS ===\n")
print(top200_overlap)

cat("\nWrote files:\n")
cat(" - comparison_run_summary.csv\n")
cat(" - comparison_top100_overlap_between_runs.csv\n")
cat(" - comparison_top200_overlap_between_runs.csv\n")
for (run_name in run_dirs) {
  cat(sprintf(" - comparison_%s_top20_final_links.csv\n", run_name))
  cat(sprintf(" - comparison_%s_high_distance_high_rank_links.csv\n", run_name))
  cat(sprintf(" - comparison_%s_mid_distance_high_rank_links.csv\n", run_name))
}