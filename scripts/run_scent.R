# ============================================================
# SCENT method
# ============================================================
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

  cand <- merge(
    gene_tss_local[, .(gene, gene_chr, tss)],
    peak_dt[, .(peak, peak_chr, peak_mid)],
    by.x = "gene_chr", by.y = "peak_chr",
    allow.cartesian = TRUE
  )
  msg("Initial cis candidate rows before distance filter: %d", nrow(cand))

  cand[, distance_bp := abs(peak_mid - tss)]
  cand <- cand[distance_bp <= opt$link_distance, .(gene, peak)]
  cand <- unique(cand)

  cand[, gene := as.character(gene)]
  cand[, peak := normalize_peak(as.character(peak))]

  msg("Cis candidate rows after distance filter: %d", nrow(cand))
  msg("Unique cis candidate genes: %d", uniqueN(cand$gene))
  msg("Unique cis candidate peaks: %d", uniqueN(cand$peak))

  candidate_universe <- unique(rbind(
    rerank_dt[, .(gene, peak)],
    linkpeaks_dt[, .(gene, peak)]
  ))
  candidate_universe[, gene := as.character(gene)]
  candidate_universe[, peak := normalize_peak(as.character(peak))]

  msg("Shared candidate-universe rows: %d", nrow(candidate_universe))
  msg("Shared candidate-universe genes: %d", uniqueN(candidate_universe$gene))
  msg("Shared candidate-universe peaks: %d", uniqueN(candidate_universe$peak))

  msg("Gene overlap before merge: %d", length(intersect(cand$gene, candidate_universe$gene)))
  msg("Peak overlap before merge: %d", length(intersect(cand$peak, candidate_universe$peak)))

  cand_shared <- merge(cand, candidate_universe, by = c("gene", "peak"))
  cand_shared <- unique(cand_shared)

  msg("Candidate rows after shared-universe intersection: %d", nrow(cand_shared))
  msg("Unique genes after shared-universe intersection: %d", uniqueN(cand_shared$gene))
  msg("Unique peaks after shared-universe intersection: %d", uniqueN(cand_shared$peak))

  if (nrow(cand_shared) == 0) {
    msg("Shared-universe intersection is empty. Falling back to native cis candidate set for SCENT.")
    cand_use <- cand
  } else {
    cand_use <- cand_shared
  }

  msg("Final SCENT candidate rows used: %d", nrow(cand_use))
  msg("Final SCENT unique genes used: %d", uniqueN(cand_use$gene))
  msg("Final SCENT unique peaks used: %d", uniqueN(cand_use$peak))

  if (nrow(cand_use) == 0) {
    stop("SCENT candidate set is empty even after fallback.")
  }

  scent_obj <- CreateSCENTObj(
    rna = rna_counts,
    atac = atac_counts,
    meta.data = as.data.frame(meta),
    peak.info = as.data.frame(cand_use[, .(gene, peak)]),
    covariates = c("log_nUMI", "log_nATAC"),
    celltypes = "celltype"
  )

  target_celltype <- if (!is.na(opt$scoring_celltype)) opt$scoring_celltype else unique(meta$celltype)[1]
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

  if ("boot_basic_p" %in% names(res)) {
    res[, score := -log10(pmax(as.numeric(boot_basic_p), 1e-300)) * sign(as.numeric(beta))]
    msg("Using SCENT score from boot_basic_p and beta.")
  } else if ("p" %in% names(res)) {
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

  msg("SCENT rows after cleaning: %d", nrow(out))
  msg("SCENT unique genes after cleaning: %d", uniqueN(out$gene))
  msg("SCENT unique peaks after cleaning: %d", uniqueN(out$peak))

  setorder(out, peak, gene, -score)
  out <- out[, .SD[1], by = .(peak, gene, method)]
  setorder(out, -score)

  msg("SCENT rows after deduplication: %d", nrow(out))
  msg("SCENT final output complete.")

  out
}
