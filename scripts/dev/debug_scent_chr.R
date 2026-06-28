suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(Signac)
  library(GenomeInfoDb)
  library(GenomicRanges)
  library(EnsDb.Hsapiens.v86)
  library(SCENT)
})

# ============================================================
# Mock opt for line-by-line debugging
# ============================================================
opt <- list(
  data_dir = "data",
  h5_file = "filtered_feature_bc_matrix.h5",
  frag_file = "atac_fragments.tsv.gz",

  output_dir = file.path("results", "debug_scent_chr22_tiny_25kb_frac020_500cells"),

  test_chr = "chr22",
  genome = "hg38",

  link_distance = 100000L,
  min_pair_frac = 0.20,
  n_cells_test = 1000L,

  scoring_celltype = NA_character_,
  scent_cores = 1L,
  scent_regr = "poisson",

  use_shared_universe = FALSE,
  seed = 42L
)

set.seed(opt$seed)
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) cat(sprintf(...), "\n")

normalize_peak <- function(x) {
  x <- as.character(x)
  x <- gsub("_", "-", x, fixed = TRUE)
  x
}

parse_peak_table <- function(peak_vec) {
  peak_vec <- normalize_peak(peak_vec)
  peak_vec2 <- gsub(":", "-", peak_vec, fixed = TRUE)
  parts <- strsplit(peak_vec2, "-", fixed = TRUE)

  ok <- lengths(parts) == 3
  if (!all(ok)) {
    bad <- unique(peak_vec[!ok])
    stop("Could not parse peak coordinates for: ", paste(head(bad, 5), collapse = ", "))
  }

  data.frame(
    peak = peak_vec,
    peak_chr = vapply(parts, `[`, character(1), 1),
    peak_start = as.numeric(vapply(parts, `[`, character(1), 2)),
    peak_end = as.numeric(vapply(parts, `[`, character(1), 3)),
    stringsAsFactors = FALSE
  )
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

  tx_df <- tx_df[order(tx_df$gene, -tx_df$is_pc, -tx_df$tx_len), ]
  tx_df <- tx_df[!duplicated(tx_df$gene), ]

  data.table(
    gene = tx_df$gene,
    gene_chr = tx_df$gene_chr,
    tss = tx_df$tss,
    tx_id = tx_df$tx_id,
    tx_biotype = tx_df$tx_biotype
  )
}


# ============================================================
# Config for chromosome sweep
# ============================================================

opt$output_dir <- file.path("results", "scent_chr_sweep_100kb_frac020_1000cells")
opt$link_distance <- 100000L
opt$min_pair_frac <- 0.20
opt$n_cells_test <- 1000L
opt$scent_cores <- 1L
opt$scent_regr <- "poisson"
opt$max_scent_candidates <- 1500L
opt$skip_existing <- TRUE

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

chroms <- c(
  "chr22", "chr21", "chr20", "chr19", "chr18", "chr17", "chr16", "chr15",
  "chr14", "chr13", "chr12", "chr11", "chr10", "chr9", "chr8", "chr7",
  "chr6", "chr5", "chr4", "chr3", "chr2", "chr1",
  "chrX", "chrY"
)
# ============================================================
# Read full data once
# ============================================================

h5_path <- file.path(opt$data_dir, opt$h5_file)
frag_path <- file.path(opt$data_dir, opt$frag_file)

stopifnot(file.exists(h5_path))
stopifnot(file.exists(frag_path))
stopifnot(file.exists(paste0(frag_path, ".tbi")))

msg("Reading 10x multiome data...")
data10x <- Read10X_h5(h5_path)

rna_counts_full <- data10x[["Gene Expression"]]
atac_counts_full <- data10x[["Peaks"]]

msg("RNA full matrix: %d genes x %d cells", nrow(rna_counts_full), ncol(rna_counts_full))
msg("ATAC full matrix: %d peaks x %d cells", nrow(atac_counts_full), ncol(atac_counts_full))

peak_info_full <- as.data.table(parse_peak_table(rownames(atac_counts_full)))
gene_tss_all <- build_gene_tss_table()

set.seed(opt$seed)
keep_cells <- sample(colnames(rna_counts_full), min(opt$n_cells_test, ncol(rna_counts_full)))

rna_counts_full <- rna_counts_full[, keep_cells, drop = FALSE]
atac_counts_full <- atac_counts_full[, keep_cells, drop = FALSE]

msg("Cell subset retained globally: %d cells", length(keep_cells))

# ============================================================
# SCENT per chromosome
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

  # Rebuild chr-specific inputs from full matrices
  chr_peak_ids <- peak_info_full[peak_chr == chr, peak]

  if (length(chr_peak_ids) == 0) {
    stop(sprintf("No ATAC peaks found for %s", chr))
  }

  rna_counts <- rna_counts_full
  atac_counts <- atac_counts_full[chr_peak_ids, , drop = FALSE]
  gene_tss <- gene_tss_all[gene_chr == chr]

  meta <- data.table(
    cell = colnames(rna_counts),
    nUMI = Matrix::colSums(rna_counts),
    nATAC = Matrix::colSums(atac_counts)
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

  # Prevalence filters
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

  # Critical fix: subset matrices to candidate genes/peaks only
  needed_genes <- intersect(unique(cand$gene), rownames(rna_counts))
  needed_peaks <- intersect(unique(cand$peak), rownames(atac_counts))

  cand <- cand[gene %in% needed_genes & peak %in% needed_peaks]

  rna_use <- rna_counts[needed_genes, , drop = FALSE]
  atac_use <- atac_counts[needed_peaks, , drop = FALSE]

  stopifnot(all(cand$gene %in% rownames(rna_use)))
  stopifnot(all(cand$peak %in% rownames(atac_use)))

  msg("Final SCENT object input check passed.")
  msg("SCENT RNA matrix after candidate subsetting: %d genes x %d cells",
      nrow(rna_use), ncol(rna_use))
  msg("SCENT ATAC matrix after candidate subsetting: %d peaks x %d cells",
      nrow(atac_use), ncol(atac_use))
  msg("SCENT candidate rows after matrix matching: %d", nrow(cand))

  if (nrow(cand) == 0) {
    stop(sprintf("SCENT candidate set is empty after matrix matching for %s", chr))
  }

  fwrite(cand, cand_csv)

  scent_obj <- CreateSCENTObj(
    rna = rna_use,
    atac = atac_use,
    meta.data = as.data.frame(meta),
    peak.info = as.data.frame(cand[, .(gene, peak)]),
    covariates = c("log_nUMI", "log_nATAC"),
    celltypes = "celltype"
  )

  target_celltype <- if (!is.na(opt$scoring_celltype)) {
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

  out <- res[, .(
    peak = normalize_peak(as.character(peak)),
    gene = as.character(gene),
    score = as.numeric(score)
  )]

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

source("scripts/run_scent_for_chr.r")

# ============================================================
# Run chromosome sweep
# ============================================================

summary_file <- file.path(opt$output_dir, "scent_chr_sweep_summary.csv")
summary_rows <- list()

for (chr in chroms) {
  res <- tryCatch(
    run_scent_for_chr(chr),
    error = function(e) {
      msg("FAILED chromosome %s: %s", chr, conditionMessage(e))

      data.table(
        chr = chr,
        status = paste0("failed: ", conditionMessage(e)),
        runtime_minutes = NA_real_,
        n_rows = NA_integer_,
        n_genes = NA_integer_,
        n_peaks = NA_integer_,
        n_positive = NA_integer_,
        n_negative = NA_integer_,
        score_min = NA_real_,
        score_median = NA_real_,
        score_max = NA_real_
      )
    }
  )

  summary_rows[[chr]] <- res

  summary_dt <- rbindlist(summary_rows, fill = TRUE)
  fwrite(summary_dt, summary_file)

  msg("Updated sweep summary: %s", summary_file)

  gc()
}

msg("SCENT chromosome sweep complete.")
msg("Summary written to: %s", summary_file)