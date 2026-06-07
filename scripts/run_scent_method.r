
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