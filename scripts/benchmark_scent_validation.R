#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(GenomicRanges)
  library(ggplot2)
})

# ============================================================
# CLI
# ============================================================
option_list <- list(
  make_option(c("--dataset"), type = "character", default = "pbmc",
              help = "Dataset/run name used in ranking filenames [default %default]"),
  make_option(c("--rankings-dir"), type = "character", default = file.path("results", "pbmc", "rankings"),
              help = "Directory containing ranking subdirectories [default %default]"),
  make_option(c("--scent-sweep-dir"), type = "character", default = file.path("results", "scent_chr_sweep_100kb_frac020_1000cells"),
              help = "Directory containing per-chromosome SCENT output CSVs [default %default]"),
  make_option(c("--output-dir"), type = "character", default = file.path("results", "pbmc", "scent_validation"),
              help = "Output directory for SCENT validation artifacts [default %default]"),
  make_option(c("--methods"), type = "character", default = "linkpeaks,coactivity,coactivity_tf,full_moddist_lambda_0_1,distance_only",
              help = "Comma-separated ranking modes to compare against SCENT [default %default]"),
  make_option(c("--top-n-values"), type = "character", default = "50,100,200",
              help = "Comma-separated top-N cutoffs [default %default]"),
  make_option(c("--top-k-compare"), type = "integer", default = 200,
              help = "Main top-K used for selected plots/tables [default %default]"),
  make_option(c("--reciprocal-overlap"), type = "double", default = 0.50,
              help = "Minimum reciprocal peak overlap for SCENT support [default %default]"),
  make_option(c("--distance-d0"), type = "double", default = 50000,
              help = "Distance scale used only if distance_score must be recomputed [default %default]"),
  make_option(c("--distal-threshold"), type = "integer", default = 50000,
              help = "Threshold in bp to define distal links [default %default]"),
  make_option(c("--restrict-to-scent-chrs"), action = "store_true", default = TRUE,
              help = "Restrict ranked methods to chromosomes present in SCENT outputs [default %default]"),
  make_option(c("--no-restrict-to-scent-chrs"), action = "store_false", dest = "restrict_to_scent_chrs",
              help = "Do not restrict ranked methods to SCENT chromosomes"),
  make_option(c("--write-standardized-links"), action = "store_true", default = FALSE,
              help = "Write full standardized per-method tables with SCENT support labels [default %default]"),
  make_option(c("--done-file"), type = "character", default = NULL,
              help = "Optional sentinel file touched on success"),
  make_option(c("--seed"), type = "integer", default = 42,
              help = "Random seed [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
names(opt) <- gsub("-", "_", names(opt), fixed = TRUE)
set.seed(opt$seed)
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

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

parse_csv_ints <- function(x) {
  vals <- as.integer(trimws(unlist(strsplit(x, ",", fixed = TRUE))))
  vals <- vals[!is.na(vals)]
  sort(unique(vals))
}

parse_csv_chars <- function(x) {
  vals <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  vals[nzchar(vals)]
}

# ============================================================
# Peak parsing / canonicalization
# ============================================================
parse_peak_table <- function(peak_vec) {
  raw_peak <- as.character(peak_vec)
  peak_tmp <- gsub("_", "-", raw_peak, fixed = TRUE)
  peak_tmp <- gsub(":", "-", peak_tmp, fixed = TRUE)
  parts <- strsplit(peak_tmp, "-", fixed = TRUE)

  ok <- lengths(parts) == 3
  if (!all(ok)) {
    bad <- unique(raw_peak[!ok])
    stop("Could not parse peak coordinates for: ", paste(head(bad, 10), collapse = ", "))
  }

  chr <- vapply(parts, `[`, character(1), 1)
  start <- suppressWarnings(as.integer(vapply(parts, `[`, character(1), 2)))
  end <- suppressWarnings(as.integer(vapply(parts, `[`, character(1), 3)))

  if (anyNA(start) || anyNA(end)) {
    bad <- raw_peak[is.na(start) | is.na(end)]
    stop("Peak start/end parse produced NA for: ", paste(head(bad, 10), collapse = ", "))
  }

  data.table(
    peak_raw = raw_peak,
    peak = paste0(chr, ":", start, "-", end),
    peak_chr = chr,
    peak_start = start,
    peak_end = end,
    peak_mid = (start + end) / 2
  )
}

canonical_peak <- function(x) {
  parse_peak_table(x)$peak
}

make_peak_granges <- function(dt) {
  p <- parse_peak_table(dt$peak)
  gr <- GRanges(seqnames = p$peak_chr, ranges = IRanges(start = p$peak_start, end = p$peak_end))
  names(gr) <- p$peak
  gr
}

add_distance_bins <- function(dt) {
  dt <- copy(dt)
  dt[, distance_bin := as.character(cut(
    as.numeric(distance_bp),
    breaks = c(-Inf, 10000, 50000, 200000, 500000, Inf),
    labels = c("0_10kb", "10_50kb", "50_200kb", "200_500kb", "gt500kb"),
    right = TRUE
  ))]
  dt[is.na(distance_bin), distance_bin := "unknown"]
  dt[]
}

# ============================================================
# Input readers
# ============================================================
read_ranked_method <- function(rankings_dir, dataset, mode) {
  path <- file.path(rankings_dir, mode, sprintf("%s_%s_ranked_links.csv", dataset, mode))
  require_file(path)

  x <- fread(path)
  require_cols(x, c("peak", "gene", "model_score"), sprintf("ranked file for %s", mode))

  out <- x[, .(
    peak = canonical_peak(peak),
    gene = as.character(gene),
    score = as.numeric(model_score),
    original_score_mode = if ("score_mode" %in% names(x)) as.character(score_mode) else mode,
    link_score = if ("link_score" %in% names(x)) as.numeric(link_score) else NA_real_,
    mul_weigh = if ("mul_weigh" %in% names(x)) as.numeric(mul_weigh) else NA_real_,
    distance_score = if ("distance_score" %in% names(x)) as.numeric(distance_score) else NA_real_,
    distance_bp = if ("distance_bp" %in% names(x)) as.numeric(distance_bp) else NA_real_,
    tf_score = if ("tf_score" %in% names(x)) as.numeric(tf_score) else NA_real_,
    distance_modifier = if ("distance_modifier" %in% names(x)) as.numeric(distance_modifier) else NA_real_,
    rank_model = if ("rank_model" %in% names(x)) as.numeric(rank_model) else NA_real_,
    rank_link = if ("rank_link" %in% names(x)) as.numeric(rank_link) else NA_real_
  )]

  out <- out[!is.na(peak) & !is.na(gene) & nzchar(peak) & nzchar(gene) & is.finite(score)]
  out[, method := mode]

  setorder(out, peak, gene, -score)
  out <- out[, .SD[1], by = .(peak, gene, method)]
  setorder(out, -score)
  out[, rank := seq_len(.N)]

  out
}

read_scent_sweep <- function(scent_sweep_dir) {
  if (!dir.exists(scent_sweep_dir)) {
    stop("SCENT sweep directory not found: ", scent_sweep_dir)
  }

  files <- list.files(
    scent_sweep_dir,
    pattern = "^scent_links_chr.*\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(files) == 0) {
    stop("No scent_links_chr*.csv files found under: ", scent_sweep_dir)
  }

  msg("Found %d SCENT chromosome files.", length(files))

  scent_list <- lapply(files, function(f) {
    x <- fread(f)
    require_cols(x, c("peak", "gene", "score"), paste("SCENT file", f))

    out <- x[, .(
      peak = canonical_peak(peak),
      gene = as.character(gene),
      score = as.numeric(score),
      source_file = f,
      source_chr = basename(dirname(f))
    )]

    out[!is.na(peak) & !is.na(gene) & nzchar(gene) & is.finite(score)]
  })

  scent <- rbindlist(scent_list, use.names = TRUE, fill = TRUE)
  scent[, method := "SCENT"]

  setorder(scent, peak, gene, -score)
  scent <- scent[, .SD[1], by = .(peak, gene, method)]
  setorder(scent, -score)
  scent[, rank := seq_len(.N)]

  scent[]
}

# ============================================================
# SCENT overlap labels
# ============================================================
add_scent_support <- function(dt, scent_dt, reciprocal_overlap) {
  dt <- copy(dt)
  dt[, scent_supported := FALSE]
  dt[, scent_best_score := NA_real_]
  dt[, scent_best_rank := NA_real_]
  dt[, scent_best_peak := NA_character_]
  dt[, scent_overlap_bp := NA_real_]
  dt[, scent_recip_overlap_query := NA_real_]
  dt[, scent_recip_overlap_scent := NA_real_]

  if (nrow(dt) == 0 || nrow(scent_dt) == 0) return(dt)

  query_gr <- make_peak_granges(dt)
  scent_gr <- make_peak_granges(scent_dt)

  hits <- findOverlaps(query_gr, scent_gr, ignore.strand = TRUE)
  if (length(hits) == 0) return(dt)

  qh <- queryHits(hits)
  sh <- subjectHits(hits)

  same_gene <- dt$gene[qh] == scent_dt$gene[sh]
  if (!any(same_gene)) return(dt)

  qh <- qh[same_gene]
  sh <- sh[same_gene]

  inter <- pintersect(query_gr[qh], scent_gr[sh])
  ov_width <- width(inter)
  recip_q <- ov_width / width(query_gr[qh])
  recip_s <- ov_width / width(scent_gr[sh])

  ok <- is.finite(recip_q) & is.finite(recip_s) &
    recip_q >= reciprocal_overlap & recip_s >= reciprocal_overlap

  if (!any(ok)) return(dt)

  support_hits <- data.table(
    row_id = qh[ok],
    scent_score = scent_dt$score[sh[ok]],
    scent_rank = scent_dt$rank[sh[ok]],
    scent_peak = scent_dt$peak[sh[ok]],
    overlap_bp = ov_width[ok],
    recip_query = recip_q[ok],
    recip_scent = recip_s[ok]
  )

  setorder(support_hits, row_id, -scent_score, scent_rank)
  best <- support_hits[, .SD[1], by = row_id]

  dt$scent_supported[best$row_id] <- TRUE
  dt$scent_best_score[best$row_id] <- best$scent_score
  dt$scent_best_rank[best$row_id] <- best$scent_rank
  dt$scent_best_peak[best$row_id] <- best$scent_peak
  dt$scent_overlap_bp[best$row_id] <- best$overlap_bp
  dt$scent_recip_overlap_query[best$row_id] <- best$recip_query
  dt$scent_recip_overlap_scent[best$row_id] <- best$recip_scent

  dt[]
}

filter_to_chrs <- function(dt, chrs) {
  pinfo <- unique(parse_peak_table(dt$peak)[, .(peak, peak_chr)])
  x <- merge(copy(dt), pinfo, by = "peak", all.x = TRUE, sort = FALSE)
  x <- x[peak_chr %in% chrs]
  x[, peak_chr := NULL]
  x[]
}

# ============================================================
# Summaries
# ============================================================
method_count_summary <- function(method_tables) {
  rbindlist(lapply(names(method_tables), function(nm) {
    x <- method_tables[[nm]]
    pinfo <- parse_peak_table(x$peak)
    data.table(
      method = nm,
      n_links = nrow(x),
      unique_genes = uniqueN(x$gene),
      unique_peaks = uniqueN(x$peak),
      n_chromosomes = uniqueN(pinfo$peak_chr),
      n_scent_supported = sum(x$scent_supported, na.rm = TRUE),
      frac_scent_supported = mean(x$scent_supported, na.rm = TRUE),
      score_min = min(x$score, na.rm = TRUE),
      score_median = median(x$score, na.rm = TRUE),
      score_max = max(x$score, na.rm = TRUE)
    )
  }), fill = TRUE)
}

make_topn_metrics <- function(dt, method_name, top_n, distal_threshold) {
  top_dt <- copy(dt[order(-score)][1:min(.N, top_n)])
  finite_dist <- top_dt$distance_bp[is.finite(top_dt$distance_bp)]

  if (nrow(top_dt) == 0) {
    return(data.table(
      method = method_name,
      top_n = top_n,
      n_links = 0,
      unique_genes = 0,
      unique_peaks = 0,
      n_scent_supported = 0,
      frac_scent_supported = NA_real_,
      median_scent_rank_supported = NA_real_,
      median_distance_bp = NA_real_,
      distal_frac = NA_real_,
      promoter_frac_10kb = NA_real_,
      median_score = NA_real_
    ))
  }

  supported <- top_dt[scent_supported == TRUE]

  data.table(
    method = method_name,
    top_n = top_n,
    n_links = nrow(top_dt),
    unique_genes = uniqueN(top_dt$gene),
    unique_peaks = uniqueN(top_dt$peak),
    n_scent_supported = nrow(supported),
    frac_scent_supported = nrow(supported) / nrow(top_dt),
    median_scent_rank_supported = if (nrow(supported) > 0) median(supported$scent_best_rank, na.rm = TRUE) else NA_real_,
    median_distance_bp = if (length(finite_dist)) median(finite_dist, na.rm = TRUE) else NA_real_,
    distal_frac = mean(top_dt$distance_bp > distal_threshold, na.rm = TRUE),
    promoter_frac_10kb = mean(top_dt$distance_bp <= 10000, na.rm = TRUE),
    median_score = median(top_dt$score, na.rm = TRUE)
  )
}

calc_distance_matched_enrichment <- function(dt, method_name, high_fraction = 0.10) {
  x <- add_distance_bins(copy(dt))
  x <- x[distance_bin != "unknown"]
  if (nrow(x) == 0) return(data.table())

  setorder(x, distance_bin, -score, peak, gene)
  x[, rank_in_bin := seq_len(.N), by = distance_bin]
  x[, n_in_bin := .N, by = distance_bin]
  x[, high_ranked := rank_in_bin <= pmax(1L, ceiling(n_in_bin * high_fraction))]

  out <- x[, {
    high <- .SD[high_ranked == TRUE]
    rest <- .SD[high_ranked == FALSE]
    a <- sum(high$scent_supported, na.rm = TRUE)
    b <- nrow(high) - a
    c <- sum(rest$scent_supported, na.rm = TRUE)
    d <- nrow(rest) - c
    odds_ratio <- ((a + 0.5) / (b + 0.5)) / ((c + 0.5) / (d + 0.5))
    data.table(
      method = method_name,
      distance_bin = distance_bin[1],
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
  }, by = distance_bin]

  out[]
}

calc_pair_overlap <- function(method_tables, ks) {
  methods <- names(method_tables)
  out <- list()
  idx <- 1L
  for (k in ks) {
    for (i in seq_along(methods)) {
      for (j in seq_along(methods)) {
        a <- unique(method_tables[[methods[i]]][order(-score)][1:min(.N, k), .(peak, gene)])
        b <- unique(method_tables[[methods[j]]][order(-score)][1:min(.N, k), .(peak, gene)])
        ov <- nrow(merge(a, b, by = c("peak", "gene")))
        out[[idx]] <- data.table(
          k = k,
          method_a = methods[i],
          method_b = methods[j],
          overlap_pairs = ov,
          frac_of_k = ov / k
        )
        idx <- idx + 1L
      }
    }
  }
  rbindlist(out)
}

write_top_supported_tables <- function(method_tables, output_dir, top_k_compare) {
  for (nm in names(method_tables)) {
    x <- copy(method_tables[[nm]][order(-score)])
    fwrite(head(x, top_k_compare), file.path(output_dir, sprintf("top%d_%s_links_with_scent_support.csv", top_k_compare, nm)))
    fwrite(head(x[scent_supported == TRUE], top_k_compare), file.path(output_dir, sprintf("top%d_%s_scent_supported_links.csv", top_k_compare, nm)))
  }
}

# ============================================================
# Main
# ============================================================
methods <- parse_csv_chars(opt$methods)
top_ns <- parse_csv_ints(opt$top_n_values)
top_ns <- sort(unique(c(top_ns, opt$top_k_compare)))

if (length(methods) == 0) stop("No methods supplied in --methods")
if (length(top_ns) == 0) stop("No top-N values supplied in --top-n-values")

msg("Reading SCENT sweep...")
scent_dt <- read_scent_sweep(opt$scent_sweep_dir)
scent_chrs <- sort(unique(parse_peak_table(scent_dt$peak)$peak_chr))
msg("SCENT chromosomes: %s", paste(scent_chrs, collapse = ", "))

msg("Reading ranked method outputs...")
method_tables <- setNames(lapply(methods, function(m) {
  read_ranked_method(opt$rankings_dir, opt$dataset, m)
}), methods)

if (isTRUE(opt$restrict_to_scent_chrs)) {
  msg("Restricting ranked methods to SCENT chromosomes.")
  method_tables <- lapply(method_tables, filter_to_chrs, chrs = scent_chrs)
}

msg("Adding SCENT overlap support labels, reciprocal overlap >= %.3f", opt$reciprocal_overlap)
method_tables <- lapply(method_tables, function(x) {
  add_scent_support(x, scent_dt, reciprocal_overlap = opt$reciprocal_overlap)
})

# Keep only finite same-chromosome distances when existing distance annotations are available.
method_tables <- lapply(method_tables, function(x) {
  if ("distance_bp" %in% names(x)) x[is.finite(distance_bp)] else x
})

combined <- rbindlist(method_tables, use.names = TRUE, fill = TRUE)
fwrite(combined, file.path(opt$output_dir, "scent_validation_all_ranked_methods_combined.csv"))

if (isTRUE(opt$write_standardized_links)) {
  for (nm in names(method_tables)) {
    fwrite(method_tables[[nm]], file.path(opt$output_dir, sprintf("scent_validation_%s_links.csv", nm)))
  }
}

counts <- method_count_summary(method_tables)
fwrite(counts, file.path(opt$output_dir, "scent_validation_method_counts.csv"))

topn_summary <- rbindlist(lapply(names(method_tables), function(nm) {
  rbindlist(lapply(top_ns, function(n) {
    make_topn_metrics(method_tables[[nm]], nm, n, opt$distal_threshold)
  }))
}), fill = TRUE)
fwrite(topn_summary, file.path(opt$output_dir, "scent_validation_topN_support_summary.csv"))

distance_enrichment <- rbindlist(lapply(names(method_tables), function(nm) {
  calc_distance_matched_enrichment(method_tables[[nm]], nm, high_fraction = 0.10)
}), fill = TRUE)
fwrite(distance_enrichment, file.path(opt$output_dir, "scent_validation_distance_matched_enrichment.csv"))

pair_overlap <- calc_pair_overlap(method_tables, ks = top_ns)
fwrite(pair_overlap, file.path(opt$output_dir, "scent_validation_pair_overlap_between_rankings.csv"))

supported_rank_summary <- rbindlist(lapply(names(method_tables), function(nm) {
  x <- method_tables[[nm]]
  supp <- x[scent_supported == TRUE]
  data.table(
    method = nm,
    n_links = nrow(x),
    n_scent_supported = nrow(supp),
    frac_scent_supported = if (nrow(x) > 0) nrow(supp) / nrow(x) else NA_real_,
    median_rank_scent_supported = if (nrow(supp) > 0) median(supp$rank, na.rm = TRUE) else NA_real_,
    mean_rank_scent_supported = if (nrow(supp) > 0) mean(supp$rank, na.rm = TRUE) else NA_real_,
    median_distance_scent_supported = if (nrow(supp) > 0) median(supp$distance_bp, na.rm = TRUE) else NA_real_
  )
}), fill = TRUE)
fwrite(supported_rank_summary, file.path(opt$output_dir, "scent_validation_supported_rank_summary.csv"))

write_top_supported_tables(method_tables, opt$output_dir, opt$top_k_compare)

manifest <- data.table(
  output = c(
    "scent_validation_method_counts.csv",
    "scent_validation_topN_support_summary.csv",
    "scent_validation_distance_matched_enrichment.csv",
    "scent_validation_supported_rank_summary.csv",
    "scent_validation_pair_overlap_between_rankings.csv",
    "scent_validation_all_ranked_methods_combined.csv"
  ),
  purpose = c(
    "Method-level size and overall SCENT support counts",
    "Top-N SCENT support, distance, distal fraction, and gene diversity summaries",
    "Within-distance-bin top-decile vs rest SCENT support enrichment",
    "Rank position summaries for SCENT-supported candidate links",
    "Exact pair overlap between ranked method top-N sets",
    "Combined long table of all ranked method candidates with SCENT support labels"
  )
)
fwrite(manifest, file.path(opt$output_dir, "scent_validation_manifest.csv"))

# ============================================================
# Plots
# ============================================================
plot_top_n <- opt$top_k_compare
plot_dt <- topn_summary[top_n == plot_top_n]

if (nrow(plot_dt) > 0) {
  p1 <- ggplot(plot_dt, aes(x = reorder(method, frac_scent_supported), y = frac_scent_supported)) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    labs(
      title = sprintf("SCENT-supported fraction among top %d links", plot_top_n),
      x = NULL,
      y = "Fraction SCENT-supported"
    )
  ggsave(file.path(opt$output_dir, "scent_validation_topK_supported_fraction.png"), p1, width = 9, height = 6, dpi = 300)

  p2 <- ggplot(plot_dt, aes(x = reorder(method, median_distance_bp), y = median_distance_bp)) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    labs(
      title = sprintf("Median distance among top %d links", plot_top_n),
      x = NULL,
      y = "Median distance to TSS (bp)"
    )
  ggsave(file.path(opt$output_dir, "scent_validation_topK_median_distance.png"), p2, width = 9, height = 6, dpi = 300)
}

if (nrow(distance_enrichment) > 0) {
  p3 <- ggplot(distance_enrichment, aes(x = distance_bin, y = method, fill = odds_ratio_high_vs_rest)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", odds_ratio_high_vs_rest)), size = 3) +
    theme_bw() +
    labs(
      title = "Distance-matched SCENT support enrichment",
      x = "Distance bin",
      y = NULL,
      fill = "OR"
    )
  ggsave(file.path(opt$output_dir, "scent_validation_distance_matched_enrichment_heatmap.png"), p3, width = 11, height = 6, dpi = 300)
}

msg("")
msg("SCENT validation complete.")
msg("Methods: %s", paste(names(method_tables), collapse = ", "))
msg("Output directory: %s", opt$output_dir)
msg("Main summary: %s", file.path(opt$output_dir, "scent_validation_topN_support_summary.csv"))
msg("Distance-matched summary: %s", file.path(opt$output_dir, "scent_validation_distance_matched_enrichment.csv"))

if (!is.null(opt$done_file)) {
  dir.create(dirname(opt$done_file), recursive = TRUE, showWarnings = FALSE)
  file.create(opt$done_file)
  msg("Touched done file: %s", opt$done_file)
}
