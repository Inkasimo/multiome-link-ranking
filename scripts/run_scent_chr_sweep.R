#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(Matrix)
  library(GenomicRanges)
  library(EnsDb.Hsapiens.v86)
  library(SCENT)
})

# ============================================================
# CLI
# ============================================================
option_list <- list(
  make_option(c("--data-dir"), type = "character", default = "data",
              help = "Directory containing the 10x multiome h5 file [default %default]"),
  make_option(c("--h5-file"), type = "character", default = NULL,
              help = "10x multiome h5 filename inside --data-dir"),
  make_option(c("--output-dir"), type = "character",
              default = file.path("results", "pbmc", "scent_chr_sweep_100kb_frac020_1000cells"),
              help = "Output directory for per-chromosome SCENT results [default %default]"),
  make_option(c("--chromosomes"), type = "character", default = "chr1",
              help = "Comma-separated chromosomes, or 'auto' for all same-chr gene/peak chromosomes [default %default]"),
  make_option(c("--link-distance"), type = "integer", default = 100000,
              help = "Maximum TSS-to-peak midpoint distance for candidate pairs [default %default]"),
  make_option(c("--min-pair-frac"), type = "double", default = 0.02,
              help = "Minimum cell fraction with nonzero gene/peak count [default %default]"),
  make_option(c("--max-cells"), type = "integer", default = 1000,
              help = "Randomly downsample cells to this count; <=0 uses all cells [default %default]"),
  make_option(c("--max-scent-candidates"), type = "integer", default = 100000,
              help = "Fail if a chromosome has more candidate pairs than this [default %default]"),
  make_option(c("--scent-cores"), type = "integer", default = 4,
              help = "Cores passed to SCENT_algorithm [default %default]"),
  make_option(c("--scent-regr"), type = "character", default = "poisson",
              help = "Regression passed to SCENT_algorithm: usually 'poisson' or 'negbin' [default %default]"),
  make_option(c("--scoring-celltype"), type = "character", default = NA_character_,
              help = "Cell type label to score. Default uses synthetic all_cells"),
  make_option(c("--skip-existing"), action = "store_true", default = TRUE,
              help = "Skip chromosomes with existing CSV and RDS [default %default]"),
  make_option(c("--no-skip-existing"), action = "store_false", dest = "skip_existing",
              help = "Do not skip existing chromosome outputs"),
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

safe_num_col <- function(dt, nm) {
  if (nm %in% names(dt)) as.numeric(dt[[nm]]) else rep(NA_real_, nrow(dt))
}

parse_csv_chars <- function(x) {
  vals <- trimws(unlist(strsplit(as.character(x), ",", fixed = TRUE)))
  vals[nzchar(vals)]
}

# ============================================================
# Peak and TSS helpers
# ============================================================
normalize_chr <- function(chr) {
  chr <- as.character(chr)
  ifelse(grepl("^chr", chr), chr, paste0("chr", chr))
}

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

  chr <- normalize_chr(vapply(parts, `[`, character(1), 1))
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

normalize_peak <- function(x) {
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
    gene_chr = normalize_chr(as.character(seqnames(tx))),
    tx_id = mcols(tx)$tx_id,
    tx_biotype = mcols(tx)$tx_biotype,
    start = start(tx),
    end = end(tx),
    strand = as.character(strand(tx)),
    stringsAsFactors = FALSE
  )

  tx_df <- tx_df[!is.na(tx_df$gene) & nzchar(tx_df$gene), ]
  tx_df <- tx_df[grepl("^chr", tx_df$gene_chr), ]

  tx_df$tss <- ifelse(tx_df$strand == "-", tx_df$end, tx_df$start)
  tx_df$tx_len <- tx_df$end - tx_df$start + 1L
  tx_df$is_pc <- tx_df$tx_biotype == "protein_coding"

  # One representative TSS per gene for tractable cis-candidate generation.
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

get_multiome_modalities <- function(h5_path) {
  requireNamespace("Seurat", quietly = TRUE)
  x <- Seurat::Read10X_h5(h5_path)

  if (!is.list(x)) {
    stop("Read10X_h5 did not return a modality list. Cannot identify RNA and ATAC matrices.")
  }

  nms <- names(x)
  rna_idx <- which(tolower(nms) %in% c("gene expression", "rna", "expression"))
  if (length(rna_idx) == 0) rna_idx <- grep("gene|rna|expression", nms, ignore.case = TRUE)

  atac_idx <- which(tolower(nms) %in% c("peaks", "atac"))
  if (length(atac_idx) == 0) atac_idx <- grep("peak|atac", nms, ignore.case = TRUE)

  if (length(rna_idx) == 0 || length(atac_idx) == 0) {
    stop("Could not identify RNA/ATAC modalities in h5. Modalities found: ", paste(nms, collapse = ", "))
  }

  rna <- x[[rna_idx[1]]]
  atac <- x[[atac_idx[1]]]

  if (!inherits(rna, "dgCMatrix")) rna <- as(rna, "dgCMatrix")
  if (!inherits(atac, "dgCMatrix")) atac <- as(atac, "dgCMatrix")

  list(rna = rna, atac = atac, modality_names = nms)
}

# ============================================================
# Load matrices
# ============================================================
if (is.null(opt$h5_file) || is.na(opt$h5_file) || !nzchar(opt$h5_file)) {
  stop("--h5-file is required")
}

h5_path <- file.path(opt$data_dir, opt$h5_file)
require_file(h5_path)

msg("Reading multiome h5: %s", h5_path)
modalities <- get_multiome_modalities(h5_path)
rna_counts_full <- modalities$rna
atac_counts_full <- modalities$atac

msg("10x modalities found: %s", paste(modalities$modality_names, collapse = ", "))
msg("Raw RNA matrix: %d genes x %d cells", nrow(rna_counts_full), ncol(rna_counts_full))
msg("Raw ATAC matrix: %d peaks x %d cells", nrow(atac_counts_full), ncol(atac_counts_full))

# Normalize peak rownames to chr:start-end.
atac_peak_norm <- normalize_peak(rownames(atac_counts_full))
if (any(duplicated(atac_peak_norm))) {
  dup <- unique(atac_peak_norm[duplicated(atac_peak_norm)])
  stop("Duplicate ATAC peak names after normalization. Example: ", paste(head(dup, 10), collapse = ", "))
}
rownames(atac_counts_full) <- atac_peak_norm

# Align cells.
common_cells <- intersect(colnames(rna_counts_full), colnames(atac_counts_full))
if (length(common_cells) == 0) stop("RNA and ATAC matrices have no shared cells.")
common_cells <- sort(common_cells)

if (is.finite(opt$max_cells) && opt$max_cells > 0 && length(common_cells) > opt$max_cells) {
  set.seed(opt$seed)
  common_cells <- sort(sample(common_cells, opt$max_cells))
  msg("Downsampled cells to max_cells=%d", opt$max_cells)
}

rna_counts_full <- rna_counts_full[, common_cells, drop = FALSE]
atac_counts_full <- atac_counts_full[, common_cells, drop = FALSE]

stopifnot(identical(colnames(rna_counts_full), colnames(atac_counts_full)))

msg("Aligned RNA matrix: %d genes x %d cells", nrow(rna_counts_full), ncol(rna_counts_full))
msg("Aligned ATAC matrix: %d peaks x %d cells", nrow(atac_counts_full), ncol(atac_counts_full))

# Library-size covariates must be global, not chromosome-local.
# Each chromosome uses the same per-cell total RNA/ATAC depth.
nUMI_full <- Matrix::colSums(rna_counts_full)
nATAC_full <- Matrix::colSums(atac_counts_full)

peak_info_full <- unique(parse_peak_table(rownames(atac_counts_full))[, .(
  peak, peak_chr, peak_start, peak_end, peak_mid
)])

msg("Building gene TSS table...")
gene_tss_all <- build_gene_tss_table()
msg("TSS genes total: %d", nrow(gene_tss_all))

available_chrs <- sort(intersect(unique(peak_info_full$peak_chr), unique(gene_tss_all$gene_chr)))
if (length(available_chrs) == 0) stop("No overlapping chromosomes between ATAC peaks and gene TSS table.")

chromosomes <- parse_csv_chars(opt$chromosomes)
if (length(chromosomes) == 1 && tolower(chromosomes) %in% c("auto", "all")) {
  chromosomes <- available_chrs
} else {
  chromosomes <- normalize_chr(chromosomes)
}
chromosomes <- unique(chromosomes)

missing_chrs <- setdiff(chromosomes, available_chrs)
if (length(missing_chrs) > 0) {
  stop("Requested chromosomes not available in both ATAC peaks and TSS table: ", paste(missing_chrs, collapse = ", "))
}

msg("SCENT chromosomes to run: %s", paste(chromosomes, collapse = ", "))
msg("SCENT parameters: link_distance=%d, min_pair_frac=%.4f, max_cells=%d, max_scent_candidates=%d, regr=%s",
    opt$link_distance, opt$min_pair_frac, opt$max_cells, opt$max_scent_candidates, opt$scent_regr)

# ============================================================
# Per-chromosome worker
# ============================================================
run_scent_for_chr <- function(chr) {
  chr_start <- Sys.time()

  msg("")
  msg("============================================================")
  msg("Starting SCENT chromosome: %s", chr)
  msg("Start time: %s", format(chr_start))
  msg("============================================================")

  chr_dir <- file.path(opt$output_dir, chr)
  dir.create(chr_dir, recursive = TRUE, showWarnings = FALSE)

  out_csv <- file.path(chr_dir, sprintf("scent_links_%s.csv", chr))
  out_rds <- file.path(chr_dir, sprintf("scent_result_%s.rds", chr))
  cand_csv <- file.path(chr_dir, sprintf("scent_candidates_%s.csv", chr))

  if (isTRUE(opt$skip_existing) && file.exists(out_csv) && file.exists(out_rds)) {
    msg("Skipping %s because output already exists.", chr)

    existing <- fread(out_csv)

    return(data.table(
      chr = chr,
      status = "skipped_existing",
      runtime_minutes = NA_real_,
      n_rows = nrow(existing),
      n_genes = uniqueN(existing$gene),
      n_peaks = uniqueN(existing$peak),
      n_positive = sum(existing$score > 0, na.rm = TRUE),
      n_negative = sum(existing$score < 0, na.rm = TRUE),
      score_min = min(existing$score, na.rm = TRUE),
      score_median = median(existing$score, na.rm = TRUE),
      score_max = max(existing$score, na.rm = TRUE)
    ))
  }

  chr_peak_ids <- peak_info_full[peak_chr == chr, peak]

  if (length(chr_peak_ids) == 0) {
    stop(sprintf("No ATAC peaks found for %s", chr))
  }

  rna_counts <- rna_counts_full
  atac_counts <- atac_counts_full[chr_peak_ids, , drop = FALSE]
  gene_tss <- gene_tss_all[gene_chr == chr]

  meta <- data.table(
    cell = colnames(rna_counts),
    nUMI = as.numeric(nUMI_full[colnames(rna_counts)]),
    nATAC = as.numeric(nATAC_full[colnames(rna_counts)])
  )

  meta[, log_nUMI := log1p(nUMI)]
  meta[, log_nATAC := log1p(nATAC)]
  meta[, percent_mito := 0]
  meta[, celltype := "all_cells"]

  stopifnot(identical(colnames(rna_counts), colnames(atac_counts)))
  stopifnot(identical(meta$cell, colnames(rna_counts)))

  msg("RNA matrix for %s: %d genes x %d cells", chr, nrow(rna_counts), ncol(rna_counts))
  msg("ATAC matrix for %s: %d peaks x %d cells", chr, nrow(atac_counts), ncol(atac_counts))
  msg("TSS genes for %s: %d", chr, nrow(gene_tss))

  expr_frac_gene <- Matrix::rowMeans(rna_counts > 0)
  expr_frac_peak <- Matrix::rowMeans(atac_counts > 0)

  keep_genes <- names(expr_frac_gene)[expr_frac_gene >= opt$min_pair_frac]
  keep_peaks <- names(expr_frac_peak)[expr_frac_peak >= opt$min_pair_frac]

  msg("Genes passing prevalence filter: %d", length(keep_genes))
  msg("Peaks passing prevalence filter: %d", length(keep_peaks))

  if (length(keep_genes) == 0 || length(keep_peaks) == 0) {
    stop(sprintf("No genes or peaks pass prevalence filter for %s", chr))
  }

  gene_tss_local <- copy(gene_tss[gene %in% keep_genes])

  peak_tbl <- parse_peak_table(keep_peaks)
  peak_tbl$peak_mid <- (peak_tbl$peak_start + peak_tbl$peak_end) / 2
  peak_dt <- as.data.table(peak_tbl)

  gene_tss_local <- gene_tss_local[gene_chr == chr]
  peak_dt <- peak_dt[peak_chr == chr]

  msg("Candidate genes on %s: %d", chr, nrow(gene_tss_local))
  msg("Candidate peaks on %s: %d", chr, nrow(peak_dt))

  if (nrow(gene_tss_local) == 0 || nrow(peak_dt) == 0) {
    stop(sprintf("No candidate genes or peaks after chr filtering for %s", chr))
  }

  cand <- merge(
    gene_tss_local[, .(gene, gene_chr, tss)],
    peak_dt[, .(peak, peak_chr, peak_mid)],
    by.x = "gene_chr",
    by.y = "peak_chr",
    allow.cartesian = TRUE
  )

  msg("Initial cis candidate rows before distance filter: %d", nrow(cand))

  cand[, distance_bp := abs(peak_mid - tss)]
  cand <- cand[distance_bp <= opt$link_distance, .(gene, peak)]
  cand <- unique(cand)

  cand[, gene := as.character(gene)]
  cand[, peak := normalize_peak(as.character(peak))]

  msg("Final SCENT candidate rows used: %d", nrow(cand))
  msg("Final SCENT unique genes used: %d", uniqueN(cand$gene))
  msg("Final SCENT unique peaks used: %d", uniqueN(cand$peak))

  if (nrow(cand) == 0) {
    stop(sprintf("SCENT candidate set is empty for %s", chr))
  }

  if (nrow(cand) > opt$max_scent_candidates) {
    stop(sprintf(
      "Too many SCENT candidates for %s: %d > max_scent_candidates=%d",
      chr, nrow(cand), opt$max_scent_candidates
    ))
  }

  needed_genes <- intersect(unique(cand$gene), rownames(rna_counts))
  needed_peaks <- intersect(unique(cand$peak), rownames(atac_counts))

  cand <- cand[gene %in% needed_genes & peak %in% needed_peaks]

  rna_use <- rna_counts[needed_genes, , drop = FALSE]
  atac_use <- atac_counts[needed_peaks, , drop = FALSE]

  stopifnot(all(cand$gene %in% rownames(rna_use)))
  stopifnot(all(cand$peak %in% rownames(atac_use)))

  msg("Final SCENT object input check passed.")
  msg("SCENT RNA matrix after candidate subsetting: %d genes x %d cells", nrow(rna_use), ncol(rna_use))
  msg("SCENT ATAC matrix after candidate subsetting: %d peaks x %d cells", nrow(atac_use), ncol(atac_use))
  msg("SCENT candidate rows after matrix matching: %d", nrow(cand))

  if (nrow(cand) == 0) {
    stop(sprintf("SCENT candidate set is empty after matrix matching for %s", chr))
  }

  fwrite(cand, cand_csv)

  meta_df <- as.data.frame(meta)
  rownames(meta_df) <- meta_df$cell

  scent_obj <- CreateSCENTObj(
    rna = rna_use,
    atac = atac_use,
    meta.data = meta_df,
    peak.info = as.data.frame(cand[, .(gene, peak)]),
    covariates = c("log_nUMI", "log_nATAC"),
    celltypes = "celltype"
  )

  target_celltype <- if (!is.na(opt$scoring_celltype) && nzchar(opt$scoring_celltype)) {
    opt$scoring_celltype
  } else {
    unique(meta$celltype)[1]
  }

  msg("SCENT target celltype: %s", target_celltype)

  scent_obj <- SCENT_algorithm(
    object = scent_obj,
    celltype = target_celltype,
    ncores = opt$scent_cores,
    regr = opt$scent_regr,
    bin = TRUE
  )

  res <- as.data.table(scent_obj@SCENT.result)

  msg("Raw SCENT result rows: %d", nrow(res))
  msg("SCENT result columns: %s", paste(names(res), collapse = ", "))

  if (!all(c("gene", "peak") %in% names(res))) {
    stop(sprintf("SCENT result is missing gene/peak columns for %s", chr))
  }

  if ("boot_basic_p" %in% names(res) && "beta" %in% names(res)) {
    res[, score := -log10(pmax(as.numeric(boot_basic_p), 1e-300)) * sign(as.numeric(beta))]
    msg("Using SCENT score from boot_basic_p and beta.")
  } else if ("p" %in% names(res) && "beta" %in% names(res)) {
    res[, score := -log10(pmax(as.numeric(p), 1e-300)) * sign(as.numeric(beta))]
    msg("Using SCENT score from p and beta.")
  } else if ("z" %in% names(res)) {
    res[, score := as.numeric(z)]
    msg("Using SCENT score from z.")
  } else {
    stop(sprintf("Could not derive a ranking score from SCENT output for %s", chr))
  }

  out <- data.table(
    peak = normalize_peak(as.character(res$peak)),
    gene = as.character(res$gene),
    score = as.numeric(res$score),
    beta = safe_num_col(res, "beta"),
    se = safe_num_col(res, "se"),
    z = safe_num_col(res, "z"),
    p = safe_num_col(res, "p"),
    boot_basic_p = safe_num_col(res, "boot_basic_p")
  )

  out[, method := "SCENT"]
  out <- out[!is.na(peak) & !is.na(gene) & is.finite(score)]

  setorder(out, peak, gene, -score)
  out <- out[, .SD[1], by = .(peak, gene, method)]
  setorder(out, -score)

  msg("SCENT rows after cleaning/deduplication: %d", nrow(out))

  fwrite(out, out_csv)

  chr_end <- Sys.time()
  runtime <- as.numeric(difftime(chr_end, chr_start, units = "mins"))

  saveRDS(
    list(
      scent_dt = out,
      raw_result = res,
      runtime_minutes = runtime,
      chr = chr,
      opt = opt
    ),
    out_rds
  )

  msg("Finished SCENT chromosome: %s", chr)
  msg("Runtime for %s: %.2f minutes", chr, runtime)

  data.table(
    chr = chr,
    status = "ok",
    runtime_minutes = runtime,
    n_rows = nrow(out),
    n_genes = uniqueN(out$gene),
    n_peaks = uniqueN(out$peak),
    n_positive = sum(out$score > 0, na.rm = TRUE),
    n_negative = sum(out$score < 0, na.rm = TRUE),
    score_min = min(out$score, na.rm = TRUE),
    score_median = median(out$score, na.rm = TRUE),
    score_max = max(out$score, na.rm = TRUE)
  )
}

# ============================================================
# Run
# ============================================================
summary_list <- lapply(chromosomes, run_scent_for_chr)
summary_dt <- rbindlist(summary_list, use.names = TRUE, fill = TRUE)

summary_csv <- file.path(opt$output_dir, "scent_chr_sweep_summary.csv")
fwrite(summary_dt, summary_csv)

# Convenience combined table for quick inspection. Validation reads per-chr files.
chr_files <- list.files(
  opt$output_dir,
  pattern = "^scent_links_chr.*\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(chr_files) > 0) {
  combined <- rbindlist(lapply(chr_files, fread), use.names = TRUE, fill = TRUE)
  fwrite(combined, file.path(opt$output_dir, "scent_links_all_chromosomes.csv"))
}

msg("")
msg("SCENT chromosome sweep complete.")
msg("Output directory: %s", opt$output_dir)
msg("Summary: %s", summary_csv)

if (!is.null(opt$done_file) && !is.na(opt$done_file) && nzchar(opt$done_file)) {
  dir.create(dirname(opt$done_file), recursive = TRUE, showWarnings = FALSE)
  file.create(opt$done_file)
  msg("Touched done file: %s", opt$done_file)
}
