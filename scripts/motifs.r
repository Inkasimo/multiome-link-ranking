#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(TFBSTools)
  library(JASPAR2022)
  library(motifmatchr)
  library(Matrix)
})

# ---------------------------
# CLI
# ---------------------------
option_list <- list(
  make_option(c("-i", "--input"),
              type = "character",
              default = "test_scores_with_distance_FINAL.csv",
              help = "Input CSV from your first script [default %default]"),
  make_option(c("-o", "--output"),
              type = "character",
              default = "test_scores_with_motifs.csv",
              help = "Output CSV with motif columns added [default %default]"),
  make_option(c("--species"),
              type = "integer",
              default = 9606,
              help = "NCBI taxonomy ID for motifs [default %default]"),
  make_option(c("--collection"),
              type = "character",
              default = "CORE",
              help = "JASPAR collection [default %default]"),
  make_option(c("--min.score"),
              type = "character",
              default = "85%",
              help = "Motif match threshold for motifmatchr [default %default]"),
  make_option(c("--motif-penalty"),
              type = "double",
              default = 0.85,
              help = "Soft multiplier used when no motif is found [default %default]"),
  make_option(c("--motif-boost"),
              type = "double",
              default = 1.00,
              help = "Soft multiplier used when motif is found [default %default]"),
  make_option(c("--top-k-names"),
              type = "integer",
              default = 3,
              help = "How many motif names to keep per peak [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---------------------------
# Helpers
# ---------------------------
parse_peak_granges <- function(peak_vec) {
  parts <- tstrsplit(peak_vec, "-", fixed = TRUE)
  if (length(parts) != 3) {
    stop("Peak format must be chr-start-end")
  }

  gr <- GRanges(
    seqnames = parts[[1]],
    ranges = IRanges(
      start = as.integer(parts[[2]]),
      end   = as.integer(parts[[3]])
    )
  )
  names(gr) <- peak_vec
  gr
}

collapse_names <- function(x, k = 3) {
  x <- unique(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return(NA_character_)
  paste(head(x, k), collapse = ";")
}

# ---------------------------
# Load results
# ---------------------------
if (!file.exists(opt$input)) {
  stop("Input file not found: ", opt$input)
}

results <- fread(opt$input)

required_cols <- c("peak", "gene", "link_score", "final_v5")
missing_cols <- setdiff(required_cols, names(results))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

cat("Loaded rows:", nrow(results), "\n")
cat("Unique peaks:", length(unique(results$peak)), "\n")

# ---------------------------
# Build unique peak set
# ---------------------------
peak_ids <- unique(results$peak)
peak_gr <- parse_peak_granges(peak_ids)

# Keep only standard seqlevels that exist in hg38 package
common_seqlevels <- intersect(seqlevels(peak_gr), seqlevels(BSgenome.Hsapiens.UCSC.hg38))
peak_gr <- keepSeqlevels(peak_gr, common_seqlevels, pruning.mode = "coarse")
peak_gr <- peak_gr[!is.na(start(peak_gr)) & !is.na(end(peak_gr)) & width(peak_gr) > 0]

cat("Peaks retained after seqlevel filtering:", length(peak_gr), "\n")

# ---------------------------
# Get JASPAR motifs
# ---------------------------
opts <- list(
  collection = opt$collection,
  tax_group = "vertebrates",
  species = opt$species
)

pfm_list <- getMatrixSet(JASPAR2022, opts)

if (length(pfm_list) == 0) {
  stop("No JASPAR motifs returned. Check species/collection settings.")
}

cat("Motifs loaded from JASPAR package:", length(pfm_list), "\n")

motif_ids <- names(pfm_list)

motif_names <- vapply(
  pfm_list,
  function(x) {
    nm <- name(x)
    if (is.null(nm) || length(nm) == 0) return(NA_character_)
    nm
  },
  character(1)
)

# ---------------------------
# Match motifs to peaks
# ---------------------------
cat("Scanning motifs across peak regions...\n")

motif_matches <- matchMotifs(
  pwms = pfm_list,
  subject = peak_gr,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  out = "matches",
  p.cutoff = NULL,
  bg = "genome",
  min.score = opt$min.score
)

motif_mat <- motifMatches(motif_matches)

# motif_mat: motifs x peaks OR peaks x motifs depending on backend
# Normalize orientation so rows = peaks, cols = motifs
if (nrow(motif_mat) == length(pfm_list) && ncol(motif_mat) == length(peak_gr)) {
  motif_mat <- t(motif_mat)
}

if (nrow(motif_mat) != length(peak_gr)) {
  stop("Unexpected motif matrix dimensions after matching.")
}

# ---------------------------
# Summarize motif support per peak
# ---------------------------
peak_has_motif <- Matrix::rowSums(motif_mat) > 0
peak_n_motifs  <- Matrix::rowSums(motif_mat)

peak_top_names <- apply(motif_mat, 1, function(hit_row) {
  hit_idx <- which(hit_row)
  collapse_names(motif_names[hit_idx], k = opt$top_k_names)
})

peak_summary <- data.table(
  peak = names(peak_gr),
  motif_present = as.integer(peak_has_motif),
  n_motifs = as.integer(peak_n_motifs),
  motif_names = unname(peak_top_names)
)

# ---------------------------
# Merge back to results
# ---------------------------
results2 <- merge(
  results,
  peak_summary,
  by = "peak",
  all.x = TRUE,
  sort = FALSE
)

results2[is.na(motif_present), motif_present := 0L]
results2[is.na(n_motifs), n_motifs := 0L]

# Simple soft motif modifier
results2[, motif_modifier := ifelse(motif_present == 1L, opt$motif_boost, opt$motif_penalty)]

# Example integrated motif-aware score
# Keep this conservative for now
results2[, final_v5_motif := final_v5 * motif_modifier]

# Rank comparison
results2[, rank_final_v5 := frank(-final_v5, ties.method = "average")]
results2[, rank_final_v5_motif := frank(-final_v5_motif, ties.method = "average")]
results2[, rank_diff_motif := rank_final_v5 - rank_final_v5_motif]

# ---------------------------
# Save
# ---------------------------
fwrite(results2, opt$output)

cat("Saved:", opt$output, "\n")
cat("Rows with motif support:", sum(results2$motif_present == 1L, na.rm = TRUE), "\n")
cat("Unique supported peaks:", length(unique(results2$peak[results2$motif_present == 1L])), "\n")

# Quick console sanity check
cat("\nTop promoted links after motif addition:\n")
print(
  results2[order(-rank_diff_motif),
           .(gene, peak, final_v5, final_v5_motif, motif_present, n_motifs, motif_names)][1:20]
)

cat("\nTop demoted links after motif addition:\n")
print(
  results2[order(rank_diff_motif),
           .(gene, peak, final_v5, final_v5_motif, motif_present, n_motifs, motif_names)][1:20]
)