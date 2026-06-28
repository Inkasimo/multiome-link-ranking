#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# =========================
# Config
# =========================
workdir <- "/mnt/g/multiome_test/results/alpha_tf_05_after_tss"

ranked_file <- file.path(workdir, "multiome_rie_ranked_links.csv")

top_n_genes <- 6
top_peaks_per_gene <- 15
distal_bp <- 50000

# Set genes you specifically want to inspect
custom_genes <- c("PAX5", "BANK1", "TCF7", "CCL5", "MS4A1", "PRF1")

# =========================
# Helpers
# =========================
require_file <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
}

check_required_cols <- function(dt, cols, label = "table") {
  missing <- setdiff(cols, names(dt))
  if (length(missing) > 0) {
    stop(sprintf(
      "Missing required columns in %s: %s",
      label,
      paste(missing, collapse = ", ")
    ))
  }
}

pick_top_peaks_per_gene <- function(dt, n = 15) {
  dt[order(rank_final_v6)][, .SD[1:min(.N, n)], by = gene]
}

pick_top_genes_by_best_link <- function(dt, n = 6) {
  dt[order(rank_final_v6), .SD[1], by = gene][order(rank_final_v6)][1:min(.N, n)]$gene
}

# Better distal selector:
# choose genes whose best distal link is high-tier, high final score, high TF support
pick_meaningful_distal_genes <- function(dt, n = 6, distal_bp = 50000) {
  distal_dt <- dt[
    distance_bp > distal_bp & tier == "High"
  ][
    order(rank_final_v6, -tf_score, -final_v6)
  ][
    , .SD[1], by = gene
  ][
    order(rank_final_v6, -tf_score, -final_v6)
  ]

  distal_dt[1:min(.N, n)]$gene
}

make_gene_link_plot <- function(dt, title_txt, out_file) {
  if (nrow(dt) == 0) {
    message("Skipping empty plot: ", out_file)
    return(invisible(NULL))
  }

  plot_dt <- copy(dt)

  if (all(c("peak_mid", "tss") %in% names(plot_dt))) {
    plot_dt[, signed_distance := peak_mid - tss]
  } else {
    plot_dt[, signed_distance := distance_bp]
  }

  label_dt <- plot_dt[order(rank_final_v6), .SD[1:min(.N, 3)], by = gene]

  p <- ggplot(plot_dt, aes(x = signed_distance, y = final_v6)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_segment(aes(xend = signed_distance, y = 0, yend = final_v6), alpha = 0.5) +
    geom_point(aes(size = mul_weigh, color = tf_score), alpha = 0.9) +
    geom_text(
      data = label_dt,
      aes(label = motif_names),
      size = 2.5,
      vjust = -0.5,
      show.legend = FALSE
    ) +
    facet_wrap(~ gene, scales = "free_x") +
    theme_bw() +
    labs(
      title = title_txt,
      x = "Distance from TSS (bp; 0 = TSS)",
      y = "final_v6",
      size = "mul_weigh",
      color = "tf_score"
    )

  ggsave(out_file, p, width = 14, height = 8, dpi = 300)
  invisible(p)
}

# =========================
# Read data
# =========================
require_file(ranked_file)
dt <- fread(ranked_file)

check_required_cols(
  dt,
  c("gene", "peak", "mul_weigh", "distance_bp", "tf_score",
    "final_v6", "rank_final_v6", "tier", "motif_names"),
  label = ranked_file
)

# =========================
# Plot 1: top interactions
# =========================
top_genes <- pick_top_genes_by_best_link(dt, n = top_n_genes)
top_dt <- dt[gene %in% top_genes]
top_dt <- pick_top_peaks_per_gene(top_dt, top_peaks_per_gene)

make_gene_link_plot(
  top_dt,
  title_txt = sprintf("Top interactions: top %d genes by final ranking", length(top_genes)),
  out_file = file.path(workdir, "multiome_rie_top_gene_link_panels.png")
)

# =========================
# Plot 2: meaningful distal interactions
# =========================
meaningful_distal_genes <- pick_meaningful_distal_genes(
  dt,
  n = top_n_genes,
  distal_bp = distal_bp
)

meaningful_distal_dt <- dt[
  gene %in% meaningful_distal_genes & distance_bp > distal_bp
]
meaningful_distal_dt <- pick_top_peaks_per_gene(meaningful_distal_dt, top_peaks_per_gene)

make_gene_link_plot(
  meaningful_distal_dt,
  title_txt = sprintf(
    "Meaningful distal interactions: High-tier links > %s kb",
    distal_bp / 1000
  ),
  out_file = file.path(workdir, "multiome_rie_meaningful_distal_gene_link_panels.png")
)

# =========================
# Plot 3: custom genes
# =========================
present_custom_genes <- intersect(custom_genes, unique(dt$gene))

custom_dt <- dt[gene %in% present_custom_genes]
custom_dt <- pick_top_peaks_per_gene(custom_dt, top_peaks_per_gene)

make_gene_link_plot(
  custom_dt,
  title_txt = sprintf(
    "Custom genes: %s",
    paste(present_custom_genes, collapse = ", ")
  ),
  out_file = file.path(workdir, "multiome_rie_custom_gene_link_panels.png")
)

# =========================
# Plot 4: custom genes, distal only
# =========================
custom_distal_dt <- dt[
  gene %in% present_custom_genes & distance_bp > distal_bp
]
custom_distal_dt <- pick_top_peaks_per_gene(custom_distal_dt, top_peaks_per_gene)

make_gene_link_plot(
  custom_distal_dt,
  title_txt = sprintf(
    "Custom genes, distal only (> %s kb): %s",
    distal_bp / 1000,
    paste(present_custom_genes, collapse = ", ")
  ),
  out_file = file.path(workdir, "multiome_rie_custom_distal_gene_link_panels.png")
)

cat("\nWrote:\n")
cat(" -", file.path(workdir, "multiome_rie_top_gene_link_panels.png"), "\n")
cat(" -", file.path(workdir, "multiome_rie_meaningful_distal_gene_link_panels.png"), "\n")
cat(" -", file.path(workdir, "multiome_rie_custom_gene_link_panels.png"), "\n")
cat(" -", file.path(workdir, "multiome_rie_custom_distal_gene_link_panels.png"), "\n")

if (length(setdiff(custom_genes, present_custom_genes)) > 0) {
  cat("\nCustom genes not found in ranked table:\n")
  cat(" -", paste(setdiff(custom_genes, present_custom_genes), collapse = ", "), "\n")
}