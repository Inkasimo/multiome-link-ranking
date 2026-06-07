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

  link_distance = 50000L,
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


####Read data

h5_path <- file.path(opt$data_dir, opt$h5_file)
frag_path <- file.path(opt$data_dir, opt$frag_file)

stopifnot(file.exists(h5_path))
stopifnot(file.exists(frag_path))
stopifnot(file.exists(paste0(frag_path, ".tbi")))

msg("Reading 10x multiome data...")
data10x <- Read10X_h5(h5_path)

rna_counts <- data10x[["Gene Expression"]]
atac_counts <- data10x[["Peaks"]]

msg("RNA matrix: %d genes x %d cells", nrow(rna_counts), ncol(rna_counts))
msg("ATAC matrix before chr filter: %d peaks x %d cells", nrow(atac_counts), ncol(atac_counts))

peak_info <- as.data.table(parse_peak_table(rownames(atac_counts)))
keep_chr_peaks <- peak_info$peak_chr == opt$test_chr

atac_counts <- atac_counts[keep_chr_peaks, , drop = FALSE]

msg(
  "ATAC matrix after %s filter: %d peaks x %d cells",
  opt$test_chr,
  nrow(atac_counts),
  ncol(atac_counts)
)

gene_tss <- build_gene_tss_table()
gene_tss <- gene_tss[gene_chr == opt$test_chr]

msg("Genes in TSS table on %s: %d", opt$test_chr, nrow(gene_tss))

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

set.seed(42)

n_cells_test <- min(opt$n_cells_test, ncol(rna_counts))
keep_cells <- sample(colnames(rna_counts), n_cells_test)

rna_counts <- rna_counts[, keep_cells, drop = FALSE]
atac_counts <- atac_counts[, keep_cells, drop = FALSE]

meta <- meta[cell %in% keep_cells]
meta <- meta[match(keep_cells, cell)]

stopifnot(identical(colnames(rna_counts), colnames(atac_counts)))
stopifnot(identical(meta$cell, colnames(rna_counts)))

msg("Cell subset retained: %d cells", ncol(rna_counts))

run_scent_method <- function() {
  msg("Running SCENT...")

  expr_frac_gene <- Matrix::rowMeans(rna_counts > 0)
  expr_frac_peak <- Matrix::rowMeans(atac_counts > 0)

  keep_genes <- names(expr_frac_gene)[expr_frac_gene >= opt$min_pair_frac]
  keep_peaks <- names(expr_frac_peak)[expr_frac_peak >= opt$min_pair_frac]

  msg("Genes passing prevalence filter: %d", length(keep_genes))
  msg("Peaks passing prevalence filter: %d", length(keep_peaks))

  gene_tss_local <- copy(gene_tss[gene %in% keep_genes])

  peak_tbl <- parse_peak_table(keep_peaks)
  peak_tbl$peak_mid <- (peak_tbl$peak_start + peak_tbl$peak_end) / 2
  peak_dt <- as.data.table(peak_tbl)

  # Extra safety: force chromosome restriction again
  gene_tss_local <- gene_tss_local[gene_chr == opt$test_chr]
  peak_dt <- peak_dt[peak_chr == opt$test_chr]

  msg("Candidate genes on %s: %d", opt$test_chr, nrow(gene_tss_local))
  msg("Candidate peaks on %s: %d", opt$test_chr, nrow(peak_dt))

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
  
  # Keep only genes/peaks actually used by SCENT candidates
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
  stop("SCENT candidate set is empty after matching to matrices.")
}

  if (nrow(cand) > 300) {
  stop(sprintf(
    "Too many SCENT candidates for overnight smoke test: %d. Reduce link_distance or increase min_pair_frac.",
    nrow(cand)
  ))
}
  fwrite(
    cand,
    file.path(opt$output_dir, sprintf("scent_candidates_%s.csv", opt$test_chr))
  )

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
    stop("SCENT result is missing gene/peak columns.")
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
    stop("Could not derive a ranking score from SCENT output.")
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

  fwrite(
    out,
    file.path(opt$output_dir, sprintf("scent_links_%s.csv", opt$test_chr))
  )

  out
}

source("scripts/run_scent_method.r")

t0 <- Sys.time()

scent_dt <- run_scent_method()

t1 <- Sys.time()
runtime <- difftime(t1, t0, units = "mins")

msg("SCENT runtime minutes: %.2f", as.numeric(runtime))

fwrite(
  scent_dt,
  file.path(opt$output_dir, sprintf("scent_links_%s_timed.csv", opt$test_chr))
)

saveRDS(
  list(
    scent_dt = scent_dt,
    runtime_minutes = as.numeric(runtime),
    opt = opt
  ),
  file.path(opt$output_dir, sprintf("scent_result_%s_timed.rds", opt$test_chr))
)