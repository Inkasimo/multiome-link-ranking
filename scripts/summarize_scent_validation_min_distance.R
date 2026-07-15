#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option(c("--input"), type = "character",
              default = "results/pbmc/scent_validation/scent_validation_all_ranked_methods_combined.csv"),
  make_option(c("--output-dir"), type = "character",
              default = "results/pbmc/scent_validation_min_distance"),
  make_option(c("--min-distances"), type = "character",
              default = "10000,25000,50000",
              help = "Comma-separated minimum distance cutoffs. Links <= cutoff are removed."),
  make_option(c("--top-n-values"), type = "character",
              default = "50,100,200,500"),
  make_option(c("--high-fraction"), type = "double",
              default = 0.10,
              help = "Top fraction within each distance bin for enrichment analysis.")
)

opt <- parse_args(OptionParser(option_list = option_list))
names(opt) <- gsub("-", "_", names(opt), fixed = TRUE)
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

parse_nums <- function(x) {
  as.numeric(trimws(unlist(strsplit(x, ",", fixed = TRUE))))
}

parse_ints <- function(x) {
  as.integer(trimws(unlist(strsplit(x, ",", fixed = TRUE))))
}

min_distances <- parse_nums(opt$min_distances)
top_ns <- parse_ints(opt$top_n_values)

x <- fread(opt$input)

required <- c("method", "peak", "gene", "score", "distance_bp", "scent_supported")
missing <- setdiff(required, names(x))
if (length(missing) > 0) {
  stop("Input missing required columns: ", paste(missing, collapse = ", "))
}

x[, score := as.numeric(score)]
x[, distance_bp := as.numeric(distance_bp)]
x[, scent_supported := tolower(as.character(scent_supported)) %in% c("true", "t", "1", "yes")]

x <- x[is.finite(score) & is.finite(distance_bp)]

add_distance_bins_fine <- function(dt) {
  dt <- copy(dt)
  dt[, distance_bin_fine := as.character(cut(
    distance_bp,
    breaks = c(-Inf, 2000, 5000, 10000, 25000, 50000, 100000, 200000, 500000, Inf),
    labels = c("0_2kb", "2_5kb", "5_10kb", "10_25kb", "25_50kb",
               "50_100kb", "100_200kb", "200_500kb", "gt500kb"),
    right = TRUE
  ))]
  dt[is.na(distance_bin_fine), distance_bin_fine := "unknown"]
  dt[]
}

method_counts <- rbindlist(lapply(min_distances, function(md) {
  y <- x[distance_bp > md]
  y[, .(
    n_links = .N,
    unique_genes = uniqueN(gene),
    unique_peaks = uniqueN(peak),
    n_scent_supported = sum(scent_supported, na.rm = TRUE),
    frac_scent_supported_all = mean(scent_supported, na.rm = TRUE),
    median_distance_bp = median(distance_bp, na.rm = TRUE)
  ), by = method][, min_distance_bp := md][]
}), fill = TRUE)

setcolorder(method_counts, c("min_distance_bp", "method"))
fwrite(method_counts, file.path(opt$output_dir, "scent_min_distance_method_counts.csv"))

topn_summary <- rbindlist(lapply(min_distances, function(md) {
  y <- x[distance_bp > md]
  rbindlist(lapply(sort(unique(y$method)), function(m) {
    ym <- y[method == m][order(-score)]
    rbindlist(lapply(top_ns, function(n) {
      top <- head(ym, n)
      data.table(
        min_distance_bp = md,
        method = m,
        top_n = n,
        n_links = nrow(top),
        unique_genes = uniqueN(top$gene),
        unique_peaks = uniqueN(top$peak),
        n_scent_supported = sum(top$scent_supported, na.rm = TRUE),
        frac_scent_supported = if (nrow(top) > 0) mean(top$scent_supported, na.rm = TRUE) else NA_real_,
        median_distance_bp = if (nrow(top) > 0) median(top$distance_bp, na.rm = TRUE) else NA_real_,
        distal_frac_50kb = if (nrow(top) > 0) mean(top$distance_bp > 50000, na.rm = TRUE) else NA_real_
      )
    }))
  }))
}), fill = TRUE)

fwrite(topn_summary, file.path(opt$output_dir, "scent_min_distance_topN_support_summary.csv"))

base <- topn_summary[method == "linkpeaks",
  .(min_distance_bp, top_n, linkpeaks_frac = frac_scent_supported)
]

delta <- merge(
  topn_summary,
  base,
  by = c("min_distance_bp", "top_n"),
  all.x = TRUE,
  sort = FALSE
)

delta[, delta_vs_linkpeaks := frac_scent_supported - linkpeaks_frac]
setorder(delta, min_distance_bp, top_n, -delta_vs_linkpeaks)

fwrite(delta, file.path(opt$output_dir, "scent_min_distance_delta_vs_linkpeaks.csv"))

distance_enrichment <- rbindlist(lapply(min_distances, function(md) {
  y <- add_distance_bins_fine(x[distance_bp > md])
  y <- y[distance_bin_fine != "unknown"]

  if (nrow(y) == 0) return(data.table())

  setorder(y, method, distance_bin_fine, -score, peak, gene)
  y[, rank_in_bin := seq_len(.N), by = .(method, distance_bin_fine)]
  y[, n_in_bin := .N, by = .(method, distance_bin_fine)]
  y[, high_ranked := rank_in_bin <= pmax(1L, ceiling(n_in_bin * opt$high_fraction))]

  y[, {
    high <- .SD[high_ranked == TRUE]
    rest <- .SD[high_ranked == FALSE]

    a <- sum(high$scent_supported, na.rm = TRUE)
    b <- nrow(high) - a
    c <- sum(rest$scent_supported, na.rm = TRUE)
    d <- nrow(rest) - c

    odds_ratio <- ((a + 0.5) / (b + 0.5)) / ((c + 0.5) / (d + 0.5))

    data.table(
      min_distance_bp = md,
      n_high = nrow(high),
      n_rest = nrow(rest),
      supported_high = a,
      supported_rest = c,
      frac_supported_high = if (nrow(high) > 0) a / nrow(high) else NA_real_,
      frac_supported_rest = if (nrow(rest) > 0) c / nrow(rest) else NA_real_,
      odds_ratio_high_vs_rest = odds_ratio,
      median_distance_high = if (nrow(high) > 0) median(high$distance_bp, na.rm = TRUE) else NA_real_,
      median_distance_rest = if (nrow(rest) > 0) median(rest$distance_bp, na.rm = TRUE) else NA_real_
    )
  }, by = .(method, distance_bin_fine)]
}), fill = TRUE)

fwrite(distance_enrichment, file.path(opt$output_dir, "scent_min_distance_distance_matched_enrichment.csv"))

cat("Wrote:\n")
cat(file.path(opt$output_dir, "scent_min_distance_method_counts.csv"), "\n")
cat(file.path(opt$output_dir, "scent_min_distance_topN_support_summary.csv"), "\n")
cat(file.path(opt$output_dir, "scent_min_distance_delta_vs_linkpeaks.csv"), "\n")
cat(file.path(opt$output_dir, "scent_min_distance_distance_matched_enrichment.csv"), "\n")
