#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
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

# Must match the SCENT sweep link_distance used to generate validation labels.
SCENT_WINDOW_BP <- 100000

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
  y <- x[distance_bp > md & distance_bp <= SCENT_WINDOW_BP]
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
  y <- x[distance_bp > md & distance_bp <= SCENT_WINDOW_BP]

  rbindlist(lapply(sort(unique(y$method)), function(m) {
    ym <- y[method == m][order(-score, peak, gene)]

    rbindlist(lapply(top_ns, function(n) {
      top <- head(ym, n)
      n_available <- nrow(ym)
      pool_fraction <- if (n_available > 0) min(n, n_available) / n_available else NA_real_

      data.table(
        min_distance_bp = md,
        method = m,
        top_n = n,
        n_available = n_available,
        pool_fraction = pool_fraction,
        pool_limited = if (!is.na(pool_fraction)) pool_fraction > 0.5 else NA,
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
  y <- add_distance_bins_fine(x[distance_bp > md & distance_bp <= SCENT_WINDOW_BP])
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

# ============================================================
# Plot: proximal-removal support fraction
# ============================================================

plot_methods <- c(
  "linkpeaks",
  "coactivity",
  "coactivity_tf",
  "full_lambda_0_1",
  "full",
  "distance_only"
)

method_labels <- c(
  "linkpeaks" = "LinkPeaks",
  "coactivity" = "Coactivity",
  "coactivity_tf" = "Coactivity + TF",
  "full_lambda_0_1" = "Full, λ = 0.1",
  "full" = "Full, λ = 0.3",
  "distance_only" = "Distance only"
)

method_colours <- c(
  "LinkPeaks" = "#4D4D4D",
  "Coactivity" = "#A6761D",
  "Coactivity + TF" = "#1B9E77",
  "Full, λ = 0.1" = "#377EB8",
  "Full, λ = 0.3" = "#984EA3",
  "Distance only" = "#E66101"
)

method_linetypes <- c(
  "LinkPeaks" = "solid",
  "Coactivity" = "longdash",
  "Coactivity + TF" = "dotdash",
  "Full, λ = 0.1" = "solid",
  "Full, λ = 0.3" = "solid",
  "Distance only" = "twodash"
)

method_shapes <- c(
  "LinkPeaks" = 16,
  "Coactivity" = 17,
  "Coactivity + TF" = 15,
  "Full, λ = 0.1" = 18,
  "Full, λ = 0.3" = 8,
  "Distance only" = 4
)

plot_dt <- copy(topn_summary[
  method %in% plot_methods &
    top_n %in% c(50L, 100L, 200L)
])

plot_dt[, min_distance_kb := min_distance_bp / 1000]

plot_dt[, method_label := method_labels[method]]
plot_dt[is.na(method_label), method_label := method]
plot_dt[, method_label := factor(method_label, levels = unname(method_labels))]

plot_dt[, top_n_label := factor(
  paste0("Top ", top_n),
  levels = paste0("Top ", sort(unique(top_n)))
)]

percent_label <- function(z) paste0(round(100 * z), "%")

p_support <- ggplot(
  plot_dt,
  aes(
    x = min_distance_kb,
    y = frac_scent_supported,
    group = method_label,
    colour = method_label,
    linetype = method_label,
    shape = method_label
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.4, stroke = 0.9) +
  facet_wrap(~ top_n_label, ncol = 3) +
  scale_x_continuous(
    breaks = c(10, 25, 50),
    labels = c("10", "25", "50"),
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  scale_y_continuous(
    labels = percent_label,
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.08))
  ) +
  scale_colour_manual(values = method_colours, drop = FALSE) +
  scale_linetype_manual(values = method_linetypes, drop = FALSE) +
  scale_shape_manual(values = method_shapes, drop = FALSE) +
  labs(
    title = "Proximal-removal stress test",
    subtitle = "SCENT support after removing candidate links within each TSS-distance threshold",
    x = "Removed links within TSS distance (kb)",
    y = "SCENT-supported fraction",
    colour = NULL,
    linetype = NULL,
    shape = NULL,
    caption = sprintf(
      "Analysis restricted to candidate links within the %d kb SCENT validation window.",
      SCENT_WINDOW_BP / 1000
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey35", fill = NA, linewidth = 0.45),
    strip.background = element_rect(fill = "grey92", color = "grey35", linewidth = 0.45),
    strip.text = element_text(face = "bold", size = 10),
    axis.text = element_text(color = "grey20"),
    axis.title = element_text(color = "grey10"),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10, margin = margin(b = 8)),
    plot.caption = element_text(size = 8.5, hjust = 0, margin = margin(t = 8)),
    legend.position = "right",
    legend.key.width = grid::unit(1.4, "lines"),
    legend.text = element_text(size = 9),
    plot.margin = margin(8, 12, 8, 8)
  )

ggsave(
  file.path(opt$output_dir, "scent_min_distance_topN_supported_fraction.png"),
  p_support,
  width = 8.5,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

cat("Wrote:\n")
cat(file.path(opt$output_dir, "scent_min_distance_method_counts.csv"), "\n")
cat(file.path(opt$output_dir, "scent_min_distance_topN_support_summary.csv"), "\n")
cat(file.path(opt$output_dir, "scent_min_distance_delta_vs_linkpeaks.csv"), "\n")
cat(file.path(opt$output_dir, "scent_min_distance_distance_matched_enrichment.csv"), "\n")
cat(file.path(opt$output_dir, "scent_min_distance_topN_supported_fraction.png"), "\n")