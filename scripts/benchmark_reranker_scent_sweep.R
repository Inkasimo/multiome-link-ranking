#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(GenomicRanges)
  library(EnsDb.Hsapiens.v86)
  library(ggplot2)
})

# ============================================================
# CLI
# ============================================================
option_list <- list(
  make_option(c("--results-dir"), type = "character",
              default = file.path("results", "alpha_tf_05_after_tss_current"),
              help = "Directory containing existing LinkPeaks/reranker outputs [default %default]"),

  make_option(c("--baseline-file"), type = "character",
              default = "multiome_rie_baseline_links_full.csv",
              help = "Baseline LinkPeaks CSV inside --results-dir [default %default]"),

  make_option(c("--ranked-file"), type = "character",
              default = "multiome_rie_ranked_links.csv",
              help = "Reranked CSV inside --results-dir [default %default]"),

  make_option(c("--test-scores-file"), type = "character",
              default = "multiome_rie_test_scores.csv",
              help = "Score table from main method inside --results-dir [default %default]"),

  make_option(c("--scent-sweep-dir"), type = "character",
              default = file.path("results", "scent_chr_sweep_100kb_frac020_1000cells"),
              help = "Directory containing per-chromosome SCENT sweep folders [default %default]"),

  make_option(c("--output-dir"), type = "character",
              default = file.path("results", "benchmark_reranker_scent_sweep"),
              help = "Output directory for benchmark artifacts [default %default]"),

  make_option(c("--distance-d0"), type = "double", default = 50000,
              help = "Distance scale for distance_score [default %default]"),

  make_option(c("--distal-threshold"), type = "integer", default = 50000,
              help = "Threshold in bp to define distal links [default %default]"),

  make_option(c("--top-k-compare"), type = "integer", default = 200,
              help = "Main top-K used for plots/overlap heatmap [default %default]"),

  make_option(c("--restrict-to-scent-chrs"), action = "store_true", default = TRUE,
              help = "Restrict all methods to chromosomes present in SCENT outputs [default %default]"),

  make_option(c("--no-restrict-to-scent-chrs"), action = "store_false",
              dest = "restrict_to_scent_chrs",
              help = "Do not restrict non-SCENT methods to SCENT chromosomes"),

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

# ============================================================
# Peak parsing / canonicalization
# ============================================================
# Canonical benchmark peak format is chr:start-end.
# This allows exact matching across files that use chr:start-end, chr-start-end, or chr_start_end.
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

build_gene_tss_table <- function() {
  tx <- ensembldb::transcripts(
    EnsDb.Hsapiens.v86,
    return.type = "GRanges",
    columns = c("gene_name", "tx_id", "tx_biotype")
  )

  tx_df <- data.frame(
    gene = mcols(tx)$gene_name,
    gene_chr = as.character(seqnames(tx)),
    tx_id = mcols(tx)$tx_id,
    tx_biotype = mcols(tx)$tx_biotype,
    start = start(tx),
    end = end(tx),
    strand = as.character(strand(tx)),
    stringsAsFactors = FALSE
  )

  tx_df <- tx_df[!is.na(tx_df$gene) & nzchar(tx_df$gene), ]
  tx_df$gene_chr <- ifelse(grepl("^chr", tx_df$gene_chr), tx_df$gene_chr, paste0("chr", tx_df$gene_chr))
  tx_df <- tx_df[grepl("^chr", tx_df$gene_chr), ]

  tx_df$tss <- ifelse(tx_df$strand == "-", tx_df$end, tx_df$start)
  tx_df$tx_len <- tx_df$end - tx_df$start + 1L
  tx_df$is_pc <- tx_df$tx_biotype == "protein_coding"

  # One TSS per gene for benchmark distance annotation.
  tx_df <- tx_df[order(tx_df$gene, -tx_df$is_pc, -tx_df$tx_len), ]
  tx_df <- tx_df[!duplicated(tx_df$gene), ]

  data.table(
    gene = as.character(tx_df$gene),
    gene_chr = as.character(tx_df$gene_chr),
    tss = as.integer(tx_df$tss),
    tx_id = as.character(tx_df$tx_id),
    tx_biotype = as.character(tx_df$tx_biotype)
  )
}

# ============================================================
# Standardization
# ============================================================
standardize_method_table <- function(dt, score_col, method_name) {
  dt <- as.data.table(dt)

  missing_cols <- setdiff(c("peak", "gene", score_col), names(dt))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "%s is missing required columns: %s",
      method_name,
      paste(missing_cols, collapse = ", ")
    ))
  }

  out <- dt[, .(
    peak = canonical_peak(peak),
    gene = as.character(gene),
    score = as.numeric(get(score_col))
  )]

  out <- out[!is.na(peak) & !is.na(gene) & nzchar(peak) & nzchar(gene) & is.finite(score)]
  out[, method := method_name]

  setorder(out, peak, gene, -score)
  out <- out[, .SD[1], by = .(peak, gene, method)]
  setorder(out, -score)

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

    missing_cols <- setdiff(c("peak", "gene", "score"), names(x))
    if (length(missing_cols) > 0) {
      stop("SCENT file missing required columns: ", f, " / ", paste(missing_cols, collapse = ", "))
    }

    chr_from_path <- basename(dirname(f))

    out <- x[, .(
      peak = canonical_peak(peak),
      gene = as.character(gene),
      score = as.numeric(score)
    )]
    out[, method := "SCENT"]
    out[, source_file := f]
    out[, source_chr := chr_from_path]

    out[!is.na(peak) & !is.na(gene) & is.finite(score)]
  })

  scent <- rbindlist(scent_list, use.names = TRUE, fill = TRUE)

  setorder(scent, peak, gene, -score)
  scent <- scent[, .SD[1], by = .(peak, gene, method)]
  setorder(scent, -score)

  scent[, c("source_file", "source_chr") := NULL]

  scent
}

add_distance_columns <- function(dt, gene_tss, distance_d0) {
  dt <- copy(as.data.table(dt))

  peak_df <- unique(parse_peak_table(dt$peak)[, .(peak, peak_chr, peak_start, peak_end, peak_mid)])
  gene_df <- gene_tss[, .(gene, gene_chr, tss)]

  dt <- merge(dt, peak_df, by = "peak", all.x = TRUE, sort = FALSE)
  dt <- merge(dt, gene_df, by = "gene", all.x = TRUE, sort = FALSE)

  dt[, distance_bp := fifelse(
    is.na(gene_chr) | is.na(peak_chr) | peak_chr != gene_chr,
    Inf,
    abs(peak_mid - tss)
  )]

  dt[, distance_score := fifelse(
    is.finite(distance_bp),
    1 / (1 + (distance_bp / distance_d0)^2),
    0
  )]

  dt[]
}

method_count_summary <- function(method_tables) {
  rbindlist(lapply(names(method_tables), function(nm) {
    x <- method_tables[[nm]]
    peak_info <- parse_peak_table(x$peak)

    data.table(
      method = nm,
      n_links = nrow(x),
      unique_genes = uniqueN(x$gene),
      unique_peaks = uniqueN(x$peak),
      n_chromosomes = uniqueN(peak_info$peak_chr),
      score_min = min(x$score, na.rm = TRUE),
      score_median = median(x$score, na.rm = TRUE),
      score_max = max(x$score, na.rm = TRUE)
    )
  }))
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
      median_distance_bp = NA_real_,
      iqr_distance_bp = NA_real_,
      distal_frac = NA_real_,
      promoter_frac_10kb = NA_real_,
      score_median = NA_real_,
      score_iqr = NA_real_
    ))
  }

  data.table(
    method = method_name,
    top_n = top_n,
    n_links = nrow(top_dt),
    unique_genes = uniqueN(top_dt$gene),
    unique_peaks = uniqueN(top_dt$peak),
    median_distance_bp = if (length(finite_dist)) median(finite_dist, na.rm = TRUE) else NA_real_,
    iqr_distance_bp = if (length(finite_dist)) IQR(finite_dist, na.rm = TRUE) else NA_real_,
    distal_frac = mean(top_dt$distance_bp > distal_threshold, na.rm = TRUE),
    promoter_frac_10kb = mean(top_dt$distance_bp <= 10000, na.rm = TRUE),
    score_median = median(top_dt$score, na.rm = TRUE),
    score_iqr = IQR(top_dt$score, na.rm = TRUE)
  )
}

calc_pair_overlap <- function(method_tables, ks = c(10, 20, 50, 100, 200, 500)) {
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

calc_gene_overlap <- function(method_tables, ks = c(10, 20, 50, 100, 200, 500)) {
  methods <- names(method_tables)
  out <- list()
  idx <- 1L

  for (k in ks) {
    for (i in seq_along(methods)) {
      for (j in seq_along(methods)) {
        a <- unique(method_tables[[methods[i]]][order(-score)][1:min(.N, k)]$gene)
        b <- unique(method_tables[[methods[j]]][order(-score)][1:min(.N, k)]$gene)
        ov <- length(intersect(a, b))

        out[[idx]] <- data.table(
          k = k,
          method_a = methods[i],
          method_b = methods[j],
          overlap_genes = ov,
          frac_of_k = ov / k
        )
        idx <- idx + 1L
      }
    }
  }

  rbindlist(out)
}

calc_spearman_shared_pairs <- function(method_tables, reference = "Reranker") {
  methods <- names(method_tables)
  out <- list()
  idx <- 1L

  if (!(reference %in% methods)) {
    stop("Reference method not found: ", reference)
  }

  ref <- method_tables[[reference]][, .(peak, gene, ref_score = score)]

  for (nm in methods) {
    if (nm == reference) next

    m <- merge(
      ref,
      method_tables[[nm]][, .(peak, gene, other_score = score)],
      by = c("peak", "gene")
    )

    out[[idx]] <- data.table(
      reference = reference,
      method = nm,
      n_shared_pairs = nrow(m),
      spearman = if (nrow(m) >= 3) {
        suppressWarnings(cor(m$ref_score, m$other_score, method = "spearman", use = "complete.obs"))
      } else {
        NA_real_
      }
    )
    idx <- idx + 1L
  }

  rbindlist(out)
}

write_top_tables <- function(method_tables, output_dir, ks = c(50, 100, 200)) {
  for (nm in names(method_tables)) {
    x <- method_tables[[nm]][order(-score)]

    for (k in ks) {
      fwrite(
        x[1:min(.N, k)],
        file.path(output_dir, sprintf("top%d_%s_links.csv", k, nm))
      )
    }
  }
}

# ============================================================
# Inputs
# ============================================================
baseline_path <- file.path(opt$results_dir, opt$baseline_file)
ranked_path <- file.path(opt$results_dir, opt$ranked_file)
test_scores_path <- file.path(opt$results_dir, opt$test_scores_file)

for (p in c(baseline_path, ranked_path, test_scores_path)) require_file(p)

msg("Reading existing reranker / baseline outputs...")
baseline_raw <- fread(baseline_path)
ranked_raw <- fread(ranked_path)
test_scores_raw <- fread(test_scores_path)

msg("Building gene TSS table...")
gene_tss <- build_gene_tss_table()

linkpeaks_dt <- standardize_method_table(baseline_raw, score_col = "score", method_name = "LinkPeaks")
rerank_dt <- standardize_method_table(ranked_raw, score_col = "final_v6", method_name = "Reranker")
coactivity_dt <- standardize_method_table(test_scores_raw, score_col = "mul_weigh", method_name = "Coactivity")
scent_dt <- read_scent_sweep(opt$scent_sweep_dir)

msg("Rows before chromosome harmonization:")
msg("  LinkPeaks:  %d", nrow(linkpeaks_dt))
msg("  Reranker:   %d", nrow(rerank_dt))
msg("  Coactivity: %d", nrow(coactivity_dt))
msg("  SCENT:      %d", nrow(scent_dt))

# Restrict all methods to chromosomes present in the SCENT sweep.
scent_peak_info <- parse_peak_table(scent_dt$peak)
scent_chrs <- sort(unique(scent_peak_info$peak_chr))

if (isTRUE(opt$restrict_to_scent_chrs)) {
  msg("Restricting all methods to SCENT chromosomes: %s", paste(scent_chrs, collapse = ", "))

  filter_to_chrs <- function(dt, chrs) {
    pinfo <- parse_peak_table(dt$peak)[, .(peak, peak_chr)]
    x <- merge(copy(dt), unique(pinfo), by = "peak", all.x = TRUE, sort = FALSE)
    x <- x[peak_chr %in% chrs]
    x[, peak_chr := NULL]
    x
  }

  linkpeaks_dt <- filter_to_chrs(linkpeaks_dt, scent_chrs)
  rerank_dt <- filter_to_chrs(rerank_dt, scent_chrs)
  coactivity_dt <- filter_to_chrs(coactivity_dt, scent_chrs)
}

msg("Rows after chromosome harmonization:")
msg("  LinkPeaks:  %d", nrow(linkpeaks_dt))
msg("  Reranker:   %d", nrow(rerank_dt))
msg("  Coactivity: %d", nrow(coactivity_dt))
msg("  SCENT:      %d", nrow(scent_dt))

# Distance-only baseline on LinkPeaks candidate universe.
distance_dt <- add_distance_columns(copy(linkpeaks_dt), gene_tss, opt$distance_d0)
distance_dt <- distance_dt[is.finite(distance_bp)]
distance_dt[, score := -distance_bp]
distance_dt[, method := "DistanceOnly"]
distance_dt <- distance_dt[, .(peak, gene, score, method)]
setorder(distance_dt, -score)

method_tables <- list(
  LinkPeaks = linkpeaks_dt,
  Reranker = rerank_dt,
  Coactivity = coactivity_dt,
  DistanceOnly = distance_dt,
  SCENT = scent_dt
)

# Add common distance columns.
method_tables <- lapply(method_tables, function(x) {
  add_distance_columns(x, gene_tss, opt$distance_d0)
})

# Remove trans/unknown distance rows from distance-dependent benchmark summaries,
# but keep all finite same-chromosome rows for fair distance metrics.
method_tables <- lapply(method_tables, function(x) x[is.finite(distance_bp)])

# Save standardized method tables.
for (nm in names(method_tables)) {
  fwrite(method_tables[[nm]], file.path(opt$output_dir, sprintf("benchmark_%s_links.csv", nm)))
}

combined <- rbindlist(method_tables, use.names = TRUE, fill = TRUE)
fwrite(combined, file.path(opt$output_dir, "benchmark_all_methods_combined.csv"))

# ============================================================
# Summaries
# ============================================================
counts <- method_count_summary(method_tables)
fwrite(counts, file.path(opt$output_dir, "benchmark_method_counts.csv"))

top_ns <- sort(unique(c(50, 100, 200, opt$top_k_compare)))
summary_all <- rbindlist(lapply(names(method_tables), function(nm) {
  rbindlist(lapply(top_ns, function(n) {
    make_topn_metrics(method_tables[[nm]], nm, n, opt$distal_threshold)
  }))
}))

fwrite(summary_all, file.path(opt$output_dir, "benchmark_summary_metrics.csv"))

pair_overlap <- calc_pair_overlap(method_tables, ks = top_ns)
gene_overlap <- calc_gene_overlap(method_tables, ks = top_ns)

fwrite(pair_overlap, file.path(opt$output_dir, "benchmark_pair_overlap.csv"))
fwrite(gene_overlap, file.path(opt$output_dir, "benchmark_gene_overlap.csv"))

spearman_vs_reranker <- calc_spearman_shared_pairs(method_tables, reference = "Reranker")
spearman_vs_linkpeaks <- calc_spearman_shared_pairs(method_tables, reference = "LinkPeaks")

fwrite(spearman_vs_reranker, file.path(opt$output_dir, "benchmark_spearman_vs_reranker.csv"))
fwrite(spearman_vs_linkpeaks, file.path(opt$output_dir, "benchmark_spearman_vs_linkpeaks.csv"))

write_top_tables(method_tables, opt$output_dir, ks = top_ns)

# Reranker-specific diagnostics.
reranker_top <- method_tables[["Reranker"]][order(-score)][1:min(.N, opt$top_k_compare)]
linkpeaks_top <- method_tables[["LinkPeaks"]][order(-score)][1:min(.N, opt$top_k_compare)]
scent_top <- method_tables[["SCENT"]][order(-score)][1:min(.N, opt$top_k_compare)]

reranker_unique_vs_linkpeaks <- fsetdiff(
  reranker_top[, .(peak, gene, score, distance_bp)],
  linkpeaks_top[, .(peak, gene, score, distance_bp)],
  all = FALSE
)

reranker_unique_vs_scent <- fsetdiff(
  reranker_top[, .(peak, gene, score, distance_bp)],
  scent_top[, .(peak, gene, score, distance_bp)],
  all = FALSE
)

fwrite(reranker_unique_vs_linkpeaks, file.path(opt$output_dir, "reranker_topK_unique_vs_linkpeaks.csv"))
fwrite(reranker_unique_vs_scent, file.path(opt$output_dir, "reranker_topK_unique_vs_scent.csv"))

# ============================================================
# Plots
# ============================================================
top_plot_n <- ifelse(opt$top_k_compare %in% summary_all$top_n, opt$top_k_compare, 200)
summary_plot <- summary_all[top_n == top_plot_n]

p1 <- ggplot(summary_plot, aes(x = reorder(method, median_distance_bp), y = median_distance_bp)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    title = sprintf("Median distance among top %d links", top_plot_n),
    x = NULL,
    y = "Median distance to TSS (bp)"
  )
ggsave(file.path(opt$output_dir, "benchmark_median_distance_topK.png"), p1, width = 9, height = 6, dpi = 300)

p2 <- ggplot(summary_plot, aes(x = reorder(method, distal_frac), y = distal_frac)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    title = sprintf("Distal fraction among top %d links (> %s bp)", top_plot_n, format(opt$distal_threshold, scientific = FALSE)),
    x = NULL,
    y = "Distal fraction"
  )
ggsave(file.path(opt$output_dir, "benchmark_distal_fraction_topK.png"), p2, width = 9, height = 6, dpi = 300)

p3 <- ggplot(summary_plot, aes(x = reorder(method, unique_genes), y = unique_genes)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    title = sprintf("Unique genes among top %d links", top_plot_n),
    x = NULL,
    y = "Unique genes"
  )
ggsave(file.path(opt$output_dir, "benchmark_unique_genes_topK.png"), p3, width = 9, height = 6, dpi = 300)

plot_dt <- rbindlist(lapply(names(method_tables), function(nm) {
  x <- copy(method_tables[[nm]][order(-score)][1:min(.N, opt$top_k_compare)])
  x[, method := nm]
  x
}), use.names = TRUE, fill = TRUE)

p4 <- ggplot(plot_dt, aes(x = distance_bp)) +
  geom_histogram(bins = 50) +
  facet_wrap(~method, scales = "free_y") +
  theme_bw() +
  labs(
    title = sprintf("Distance distribution of top %d links", opt$top_k_compare),
    x = "Distance to TSS (bp)",
    y = "Count"
  )
ggsave(file.path(opt$output_dir, "benchmark_distance_distribution_topK.png"), p4, width = 14, height = 10, dpi = 300)

ovk <- pair_overlap[k == top_plot_n]
p5 <- ggplot(ovk, aes(x = method_a, y = method_b, fill = overlap_pairs)) +
  geom_tile() +
  geom_text(aes(label = overlap_pairs), size = 3) +
  theme_bw() +
  coord_equal() +
  labs(
    title = sprintf("Top-%d exact pair overlap", top_plot_n),
    x = NULL,
    y = NULL
  )
ggsave(file.path(opt$output_dir, "benchmark_pair_overlap_heatmap.png"), p5, width = 8, height = 7, dpi = 300)

gene_ovk <- gene_overlap[k == top_plot_n]
p6 <- ggplot(gene_ovk, aes(x = method_a, y = method_b, fill = overlap_genes)) +
  geom_tile() +
  geom_text(aes(label = overlap_genes), size = 3) +
  theme_bw() +
  coord_equal() +
  labs(
    title = sprintf("Top-%d gene overlap", top_plot_n),
    x = NULL,
    y = NULL
  )
ggsave(file.path(opt$output_dir, "benchmark_gene_overlap_heatmap.png"), p6, width = 8, height = 7, dpi = 300)

# Reranker vs LinkPeaks score scatter on shared exact pairs.
lr_shared <- merge(
  method_tables[["LinkPeaks"]][, .(peak, gene, linkpeaks = score)],
  method_tables[["Reranker"]][, .(peak, gene, reranker = score, distance_bp)],
  by = c("peak", "gene")
)

if (nrow(lr_shared) > 0) {
  p7 <- ggplot(lr_shared, aes(x = linkpeaks, y = reranker, color = distance_bp)) +
    geom_point(alpha = 0.5, size = 1) +
    theme_bw() +
    labs(
      title = "Shared-pair score comparison: LinkPeaks vs Reranker",
      x = "LinkPeaks score",
      y = "Reranker score"
    )
  ggsave(file.path(opt$output_dir, "benchmark_linkpeaks_vs_reranker_scatter.png"), p7, width = 8, height = 6, dpi = 300)
}

# ============================================================
# Console summary
# ============================================================
msg("")
msg("Benchmark complete.")
msg("Compared methods: %s", paste(names(method_tables), collapse = ", "))
msg("SCENT chromosomes used: %s", paste(scent_chrs, collapse = ", "))
msg("Main summary: %s", file.path(opt$output_dir, "benchmark_summary_metrics.csv"))
msg("Method counts: %s", file.path(opt$output_dir, "benchmark_method_counts.csv"))
msg("Pair overlap: %s", file.path(opt$output_dir, "benchmark_pair_overlap.csv"))
msg("Gene overlap: %s", file.path(opt$output_dir, "benchmark_gene_overlap.csv"))
msg("Spearman vs Reranker: %s", file.path(opt$output_dir, "benchmark_spearman_vs_reranker.csv"))
msg("Standardized link tables written to: %s", opt$output_dir)
