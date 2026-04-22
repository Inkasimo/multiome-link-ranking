#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(Signac)
  library(GenomeInfoDb)
  library(GenomicRanges)
  library(SummarizedExperiment)
  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(ArchR)
  library(SCENT)
  library(GenomicRanges)
})

# ============================================================
# CLI
# ============================================================
option_list <- list(
  make_option(c("--results-dir"), type = "character", default = file.path("results", "alpha_tf_05_after_tss_current"),
              help = "Directory containing existing LinkPeaks/reranker outputs [default %default]"),
  make_option(c("--baseline-file"), type = "character", default = "multiome_rie_baseline_links_full.csv",
              help = "Baseline LinkPeaks CSV inside --results-dir [default %default]"),
  make_option(c("--ranked-file"), type = "character", default = "multiome_rie_ranked_links.csv",
              help = "Reranked method CSV inside --results-dir [default %default]"),
  make_option(c("--test-scores-file"), type = "character", default = "multiome_rie_test_scores.csv",
              help = "Score table from main method inside --results-dir [default %default]"),
  make_option(c("--data-dir"), type = "character", default = "data",
              help = "Directory containing 10x multiome input files [default %default]"),
  make_option(c("--h5-file"), type = "character", default = "filtered_feature_bc_matrix.h5",
              help = "10x h5 matrix filename inside --data-dir [default %default]"),
  make_option(c("--frag-file"), type = "character", default = "atac_fragments.tsv.gz",
              help = "ATAC fragments filename inside --data-dir [default %default]"),
  make_option(c("--output-dir"), type = "character", default = file.path("results", "benchmark_panel"),
              help = "Output directory for benchmark artifacts [default %default]"),
  make_option(c("--sample-name"), type = "character", default = "multiome_sample",
              help = "Sample name for ArchR import [default %default]"),
  make_option(c("--genome"), type = "character", default = "hg38",
              help = "Genome build for ArchR [default %default]"),
  make_option(c("--link-distance"), type = "integer", default = 500000,
              help = "Max cis distance for candidate generation / summaries [default %default]"),
  make_option(c("--distance-d0"), type = "double", default = 50000,
              help = "Distance scale for distance summaries [default %default]"),
  make_option(c("--distal-threshold"), type = "integer", default = 50000,
              help = "Threshold in bp to define distal links [default %default]"),
  make_option(c("--top-n"), type = "integer", default = 100,
              help = "Top N genes for ORA and top-N summaries [default %default]"),
  make_option(c("--top-k-compare"), type = "integer", default = 200,
              help = "Top K links for overlap/plot diagnostics [default %default]"),
  make_option(c("--min-pair-frac"), type = "double", default = 0.05,
              help = "Minimum fraction of cells for gene/peak prevalence filter in SCENT candidate list [default %default]"),
  make_option(c("--scoring-celltype"), type = "character", default = NA,
              help = "Optional celltype value for SCENT. If NA, uses all cells as one group [default %default]"),
  make_option(c("--archr-threads"), type = "integer", default = 4,
              help = "Threads for ArchR [default %default]"),
  make_option(c("--scent-cores"), type = "integer", default = 4,
              help = "Cores for SCENT [default %default]"),
  make_option(c("--scent-regr"), type = "character", default = "poisson",
              help = "Regression family for SCENT: poisson or negbin [default %default]"),
  make_option(c("--ora-show-category"), type = "integer", default = 15,
              help = "Number of categories shown in ORA dotplots [default %default]"),
  make_option(c("--seed"), type = "integer", default = 42,
              help = "Random seed [default %default]"),
  make_option(c("--skip-archr"), action = "store_true", default = FALSE,
              help = "Skip ArchR run [default %default]"),
  make_option(c("--skip-scent"), action = "store_true", default = FALSE,
              help = "Skip SCENT run [default %default]"),
  make_option(c("--archr-force"), action = "store_true", default = FALSE,
              help = "Force recreation of ArchR project files [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Helpers
# ============================================================
msg <- function(...) cat(sprintf(...), "\n")

require_file <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
}

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

safe_bitr <- function(genes) {
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (length(genes) == 0) {
    return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
  }
  out <- suppressWarnings(
    clusterProfiler::bitr(
      genes,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
  )
  unique(out)
}

safe_enrich_go <- function(entrez_ids, universe_ids) {
  if (length(entrez_ids) == 0 || length(universe_ids) == 0) return(NULL)
  tryCatch({
    enrichGO(
      gene = unique(entrez_ids),
      universe = unique(universe_ids),
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.20,
      readable = TRUE
    )
  }, error = function(e) NULL)
}

empty_plot <- function(title_txt) {
  ggplot() +
    annotate("text", x = 1, y = 1, label = "No significant enrichment", size = 5) +
    xlim(0, 2) + ylim(0, 2) +
    theme_void() +
    ggtitle(title_txt)
}

make_dotplot_safe <- function(ora_obj, title_txt, show_n = 15) {
  if (is.null(ora_obj) || nrow(as.data.frame(ora_obj)) == 0) {
    return(empty_plot(title_txt))
  }
  dotplot(ora_obj, showCategory = show_n) + ggtitle(title_txt) + theme_bw()
}

standardize_method_table <- function(dt, score_col, method_name) {
  dt <- as.data.table(dt)
  stopifnot(all(c("peak", "gene", score_col) %in% names(dt)))
  out <- dt[, .(
    peak = normalize_peak(as.character(peak)),
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

add_distance_columns <- function(dt, gene_tss) {
  dt <- copy(as.data.table(dt))
  peak_df <- unique(parse_peak_table(dt$peak))
  peak_df$peak_mid <- (peak_df$peak_start + peak_df$peak_end) / 2
  dt <- merge(dt, peak_df, by = "peak", all.x = TRUE, sort = FALSE)
  dt <- merge(dt, gene_tss[, .(gene, gene_chr, tss)], by = "gene", all.x = TRUE, sort = FALSE)
  dt[, distance_bp := ifelse(is.na(gene_chr) | peak_chr != gene_chr, Inf, abs(peak_mid - tss))]
  dt[, distance_score := ifelse(is.finite(distance_bp), 1 / (1 + (distance_bp / opt$distance_d0)^2), 0)]
  dt
}

safe_extract_counts <- function(seurat_obj, assay = "RNA", layer = "counts") {
  tryCatch({
    GetAssayData(seurat_obj, assay = assay, layer = layer)
  }, error = function(e1) {
    tryCatch({
      GetAssayData(seurat_obj, assay = assay, slot = layer)
    }, error = function(e2) {
      stop("Could not extract assay data for assay=", assay, ", layer=", layer)
    })
  })
}

make_topn_metrics <- function(dt, method_name, top_n, distal_threshold) {
  top_dt <- copy(dt[order(-score)][1:min(.N, top_n)])
  if (nrow(top_dt) == 0) {
    return(data.table(
      method = method_name,
      top_n = top_n,
      n_links = 0,
      unique_genes = 0,
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
    median_distance_bp = median(top_dt$distance_bp[is.finite(top_dt$distance_bp)], na.rm = TRUE),
    iqr_distance_bp = IQR(top_dt$distance_bp[is.finite(top_dt$distance_bp)], na.rm = TRUE),
    distal_frac = mean(top_dt$distance_bp > distal_threshold, na.rm = TRUE),
    promoter_frac_10kb = mean(top_dt$distance_bp <= 10000, na.rm = TRUE),
    score_median = median(top_dt$score, na.rm = TRUE),
    score_iqr = IQR(top_dt$score, na.rm = TRUE)
  )
}

make_ora_result <- function(dt, method_name, background_genes, top_n) {
  top_genes <- unique(dt[order(-score)][1:min(.N, top_n)]$gene)
  bg_map <- safe_bitr(background_genes)
  top_map <- safe_bitr(top_genes)
  bg_ids <- unique(bg_map$ENTREZID)
  top_ids <- unique(top_map$ENTREZID)
  ora <- safe_enrich_go(top_ids, bg_ids)
  ora_df <- if (is.null(ora)) data.frame() else as.data.frame(ora)
  list(
    method = method_name,
    genes = top_genes,
    ora = ora,
    ora_df = ora_df,
    n_terms = nrow(ora_df)
  )
}

calc_pair_overlap <- function(method_tables, ks = c(10, 20, 50, 100, 200)) {
  methods <- names(method_tables)
  out <- list()
  idx <- 1L
  for (k in ks) {
    for (i in seq_along(methods)) {
      for (j in seq_along(methods)) {
        a <- unique(method_tables[[methods[i]]][order(-score)][1:min(.N, k), .(peak, gene)])
        b <- unique(method_tables[[methods[j]]][order(-score)][1:min(.N, k), .(peak, gene)])
        ov <- nrow(merge(a, b, by = c("peak", "gene")))
        out[[idx]] <- data.table(k = k, method_a = methods[i], method_b = methods[j], overlap_pairs = ov)
        idx <- idx + 1L
      }
    }
  }
  rbindlist(out)
}

calc_gene_overlap <- function(method_tables, ks = c(10, 20, 50, 100, 200)) {
  methods <- names(method_tables)
  out <- list()
  idx <- 1L
  for (k in ks) {
    for (i in seq_along(methods)) {
      for (j in seq_along(methods)) {
        a <- unique(method_tables[[methods[i]]][order(-score)][1:min(.N, k)]$gene)
        b <- unique(method_tables[[methods[j]]][order(-score)][1:min(.N, k)]$gene)
        ov <- length(intersect(a, b))
        out[[idx]] <- data.table(k = k, method_a = methods[i], method_b = methods[j], overlap_genes = ov)
        idx <- idx + 1L
      }
    }
  }
  rbindlist(out)
}


make_peak_gr <- function(peaks) {
  dt <- as.data.table(parse_peak_table(peaks))
  GRanges(
    seqnames = dt$peak_chr,
    ranges = IRanges(dt$peak_start, dt$peak_end),
    peak = dt$peak
  )
}

map_peaks_reciprocal <- function(query_peaks, subject_peaks, min_recip = 0.5) {
  qgr <- make_peak_gr(unique(query_peaks))
  sgr <- make_peak_gr(unique(subject_peaks))

  hits <- findOverlaps(qgr, sgr, ignore.strand = TRUE)
  if (length(hits) == 0) return(data.table())

  qh <- queryHits(hits)
  sh <- subjectHits(hits)

  ov <- pintersect(qgr[qh], sgr[sh])
  ov_width <- width(ov)

  q_width <- width(qgr[qh])
  s_width <- width(sgr[sh])

  dt <- data.table(
    query_peak = mcols(qgr[qh])$peak,
    subject_peak = mcols(sgr[sh])$peak,
    overlap_bp = ov_width,
    frac_query = ov_width / q_width,
    frac_subject = ov_width / s_width
  )

  dt <- dt[
    frac_query >= min_recip & frac_subject >= min_recip
  ]

  # if multiple matches remain, keep the best one
  setorder(dt, query_peak, -pmin(frac_query, frac_subject), -overlap_bp)
  dt <- dt[, .SD[1], by = query_peak]

  dt
}
# ============================================================
# Inputs
# ============================================================
baseline_path <- file.path(opt$results_dir, opt$baseline_file)
ranked_path <- file.path(opt$results_dir, opt$ranked_file)
test_scores_path <- file.path(opt$results_dir, opt$test_scores_file)
h5_path <- file.path(opt$data_dir, opt$h5_file)
frag_path <- file.path(opt$data_dir, opt$frag_file)
frag_index_path <- paste0(frag_path, ".tbi")

for (p in c(baseline_path, ranked_path, test_scores_path, h5_path, frag_path, frag_index_path)) require_file(p)

# ============================================================
# Read existing method outputs
# ============================================================
msg("Reading existing outputs...")
baseline_raw <- fread(baseline_path)
ranked_raw <- fread(ranked_path)
test_scores_raw <- fread(test_scores_path)

gene_tss <- build_gene_tss_table()

linkpeaks_dt <- standardize_method_table(baseline_raw, score_col = "score", method_name = "LinkPeaks")
rerank_dt <- standardize_method_table(ranked_raw, score_col = "final_v6", method_name = "Reranker")
coactivity_dt <- standardize_method_table(test_scores_raw, score_col = "mul_weigh", method_name = "Coactivity")

distance_dt <- add_distance_columns(copy(linkpeaks_dt), gene_tss)
distance_dt[, score := -distance_bp]
distance_dt[, method := "DistanceOnly"]
setorder(distance_dt, -score)
distance_dt <- distance_dt[, .SD[1], by = .(peak, gene, method)]

# ============================================================
# Read raw input once for ArchR and SCENT
# ============================================================
msg("Reading 10x multiome data...")
data10x <- Read10X_h5(h5_path)
if (!all(c("Gene Expression", "Peaks") %in% names(data10x))) {
  stop("Expected Gene Expression and Peaks matrices in 10x H5.")
}
rna_counts <- data10x[["Gene Expression"]]
atac_counts <- data10x[["Peaks"]]

msg("Creating lightweight Seurat object for metadata / counts alignment...")
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- opt$genome
obj <- CreateSeuratObject(counts = rna_counts, assay = "RNA")
obj[["ATAC"]] <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = opt$genome,
  fragments = frag_path,
  annotation = annotations
)

# Basic metadata for SCENT
meta <- data.table(
  cell = colnames(obj),
  nUMI = Matrix::colSums(rna_counts),
  nATAC = Matrix::colSums(atac_counts)
)
meta[, log_nUMI := log1p(nUMI)]
meta[, log_nATAC := log1p(nATAC)]
meta[, percent_mito := 0]

if (!is.na(opt$scoring_celltype)) {
  if (!("seurat_clusters" %in% colnames(obj@meta.data))) {
    DefaultAssay(obj) <- "RNA"
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, verbose = FALSE)
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(obj, verbose = FALSE)
    DefaultAssay(obj) <- "ATAC"
    obj <- RunTFIDF(obj)
    obj <- FindTopFeatures(obj, min.cutoff = "q0")
    obj <- RunSVD(obj)
    obj <- FindMultiModalNeighbors(obj,
      reduction.list = list("pca", "lsi"),
      dims.list = list(1:30, 2:30)
    )
    obj <- FindClusters(obj, graph.name = "wsnn", resolution = 0.5, verbose = FALSE)
  }
  meta[, celltype := as.character(obj$seurat_clusters)]
} else {
  meta[, celltype := "all_cells"]
}

# ============================================================
# ArchR method
# ============================================================
run_archr_method <- function() {
  msg("Running ArchR Peak2GeneLinks...")
  archr_dir <- file.path(opt$output_dir, "archr_project")
  dir.create(archr_dir, recursive = TRUE, showWarnings = FALSE)

  addArchRThreads(threads = opt$archr_threads)
  addArchRGenome(opt$genome)
  addArchRLocking(locking = FALSE)

  arrow_files <- file.path(archr_dir, paste0(opt$sample_name, ".arrow"))
  if (!file.exists(arrow_files) || isTRUE(opt$archr_force)) {
    arrow_files <- createArrowFiles(
      inputFiles = setNames(frag_path, opt$sample_name),
      sampleNames = opt$sample_name,
      filterTSS = 0,
      filterFrags = 0,
      addTileMat = TRUE,
      addGeneScoreMat = TRUE,
      force = TRUE
    )
  }

  proj <- ArchRProject(
    ArrowFiles = arrow_files,
    outputDirectory = archr_dir,
    copyArrows = FALSE
  )

  seRNA <- import10xFeatureMatrix(
    input = setNames(h5_path, opt$sample_name),
    names = opt$sample_name
  )

  proj <- addGeneExpressionMatrix(
    input = proj,
    seRNA = seRNA,
    force = TRUE
  )

  proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    iterations = 1,
    clusterParams = list(
      resolution = 0.2,
      sampleCells = min(5000, nCells(proj)),
      n.start = 10
    ),
    varFeatures = 10000,
    dimsToUse = 1:30,
    sampleCellsPre = min(5000, nCells(proj)),
    projectCellsPre = FALSE,
    saveIterations = FALSE,
    force = TRUE
  )

  proj <- addClusters(
    input = proj,
    reducedDims = "IterativeLSI",
    name = "Clusters",
    method = "Seurat",
    sampleCells = 10000,
    force = TRUE
  )

  proj <- addGroupCoverages(
    ArchRProj = proj,
    groupBy = "Clusters",
    force = TRUE
  )

  proj <- addReproduciblePeakSet(
    ArchRProj = proj,
    groupBy = "Clusters",
    peakMethod = "Tiles",
    force = TRUE
  )

  proj <- addPeakMatrix(
    ArchRProj = proj,
    force = TRUE
  )

  proj <- addPeak2GeneLinks(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    useMatrix = "GeneExpressionMatrix",
    corCutOff = 0.75,
    maxDist = opt$link_distance
  )

  p2g <- getPeak2GeneLinks(
    ArchRProj = proj,
    corCutOff = -1,
    FDRCutOff = 1,
    varCutOffATAC = 0,
    varCutOffRNA = 0,
    returnLoops = FALSE
  )

  
  print(class(p2g))
  print(length(p2g))
  print(head(as.data.frame(p2g)))


  p2g_loops <- getPeak2GeneLinks(
    ArchRProj = proj,
    corCutOff = -1,
    FDRCutOff = 1,
    varCutOffATAC = 0,
    varCutOffRNA = 0,
    returnLoops = TRUE
  )

  print(class(p2g_loops))
  print(length(p2g_loops))
  print(head(as.data.frame(p2g_loops)))

  p2g_df <- as.data.table(as.data.frame(p2g))

  if (nrow(p2g_df) == 0) {
    warning("ArchR returned zero peak-gene links after relaxed retrieval.")
    return(data.table(
      peak = character(),
      gene = character(),
      score = numeric(),
      method = character()
    ))
  }

  peak_set <- S4Vectors::metadata(p2g)$peakSet
  gene_set <- S4Vectors::metadata(p2g)$geneSet

  peak_df <- data.table(
    idxATAC = mcols(peak_set)$idx,
    peak = paste0(
      as.character(seqnames(peak_set)),
      ":",
      start(peak_set),
      "-",
      end(peak_set)
    )
  )

  gene_df <- data.table(
    idxRNA = mcols(gene_set)$idx,
    gene = mcols(gene_set)$name
  )

  out <- merge(p2g_df, peak_df, by = "idxATAC", all.x = TRUE)
  out <- merge(out, gene_df, by = "idxRNA", all.x = TRUE)

  out <- out[, .(
    peak = normalize_peak(as.character(peak)),
    gene = as.character(gene),
    score = as.numeric(Correlation),
    method = "ArchR"
  )]

  out <- out[!is.na(peak) & !is.na(gene) & is.finite(score)]
  setorder(out, peak, gene, -score)
  out <- out[, .SD[1], by = .(peak, gene, method)]
  setorder(out, -score)
  out
}

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



archr_dt <- NULL
if (!isTRUE(opt$skip_archr)) {
  archr_dt <- tryCatch(run_archr_method(), error = function(e) {
    msg("ArchR failed: %s", conditionMessage(e))
    NULL
  })
}

scent_dt <- NULL
if (!isTRUE(opt$skip_scent)) {
  scent_dt <- tryCatch(run_scent_method(), error = function(e) {
    msg("SCENT failed: %s", conditionMessage(e))
    NULL
  })
}

# ============================================================
# Collect methods and add common annotations
# ============================================================
method_tables <- list(
  LinkPeaks = linkpeaks_dt,
  Reranker = rerank_dt,
  Coactivity = coactivity_dt,
  DistanceOnly = distance_dt[, .(peak, gene, score, method)]
)
if (!is.null(archr_dt)) method_tables[["ArchR"]] <- archr_dt
if (!is.null(scent_dt)) method_tables[["SCENT"]] <- scent_dt

method_tables <- lapply(method_tables, function(x) add_distance_columns(x, gene_tss))
combined <- rbindlist(method_tables, use.names = TRUE, fill = TRUE)

for (nm in names(method_tables)) {
  fwrite(method_tables[[nm]], file.path(opt$output_dir, sprintf("benchmark_%s_links.csv", nm)))
}

# ============================================================
# Shared background and ORA
# ============================================================
background_genes <- unique(combined$gene)
ora_results <- lapply(names(method_tables), function(nm) {
  make_ora_result(method_tables[[nm]], nm, background_genes, opt$top_n)
})
names(ora_results) <- names(method_tables)

ora_summary <- rbindlist(lapply(ora_results, function(x) {
  data.table(method = x$method, n_ora_terms = x$n_terms)
}))
fwrite(ora_summary, file.path(opt$output_dir, "benchmark_ora_summary.csv"))

for (nm in names(ora_results)) {
  fwrite(ora_results[[nm]]$ora_df, file.path(opt$output_dir, sprintf("benchmark_%s_ora_GO_BP.csv", nm)))
  p <- make_dotplot_safe(ora_results[[nm]]$ora, sprintf("ORA: %s", nm), opt$ora_show_category)
  ggsave(file.path(opt$output_dir, sprintf("benchmark_%s_ora_dotplot.png", nm)), p, width = 10, height = 7, dpi = 300)
}

# ============================================================
# Robust method summaries
# ============================================================
summary_top50 <- rbindlist(lapply(names(method_tables), function(nm) {
  make_topn_metrics(method_tables[[nm]], nm, 50, opt$distal_threshold)
}))
summary_top100 <- rbindlist(lapply(names(method_tables), function(nm) {
  make_topn_metrics(method_tables[[nm]], nm, 100, opt$distal_threshold)
}))
summary_top200 <- rbindlist(lapply(names(method_tables), function(nm) {
  make_topn_metrics(method_tables[[nm]], nm, 200, opt$distal_threshold)
}))

summary_all <- rbindlist(list(summary_top50, summary_top100, summary_top200))
summary_all <- merge(summary_all, ora_summary, by = "method", all.x = TRUE)
fwrite(summary_all, file.path(opt$output_dir, "benchmark_summary_metrics.csv"))

# Spearman against LinkPeaks and Reranker within shared pairs
cor_rows <- list()
idx <- 1L
for (nm in names(method_tables)) {
  if (nm == "LinkPeaks") next
  m <- merge(
    method_tables[["LinkPeaks"]][, .(peak, gene, score_linkpeaks = score)],
    method_tables[[nm]][, .(peak, gene, score_other = score)],
    by = c("peak", "gene")
  )
  cor_rows[[idx]] <- data.table(
    method = nm,
    n_shared_pairs = nrow(m),
    spearman_vs_linkpeaks = suppressWarnings(cor(m$score_linkpeaks, m$score_other, method = "spearman", use = "complete.obs"))
  )
  idx <- idx + 1L
}
cor_dt <- rbindlist(cor_rows)
fwrite(cor_dt, file.path(opt$output_dir, "benchmark_spearman_vs_linkpeaks.csv"))

pair_overlap <- calc_pair_overlap(method_tables, ks = c(10, 20, 50, 100, opt$top_k_compare))
gene_overlap <- calc_gene_overlap(method_tables, ks = c(10, 20, 50, 100, opt$top_k_compare))
fwrite(pair_overlap, file.path(opt$output_dir, "benchmark_pair_overlap.csv"))
fwrite(gene_overlap, file.path(opt$output_dir, "benchmark_gene_overlap.csv"))

# ============================================================
# Plots
# ============================================================
# 1. Method summary bars
p1 <- ggplot(summary_all[top_n == 50], aes(x = reorder(method, n_ora_terms), y = n_ora_terms, fill = method)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  theme_bw() +
  labs(title = "GO BP ORA terms by method (top 50 genes)", x = NULL, y = "Significant GO BP terms")

ggsave(file.path(opt$output_dir, "benchmark_ora_term_counts_top50.png"), p1, width = 9, height = 6, dpi = 300)

p2 <- ggplot(summary_all[top_n == 50], aes(x = reorder(method, median_distance_bp), y = median_distance_bp, fill = method)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  theme_bw() +
  labs(title = "Median distance among top 50 links", x = NULL, y = "Median distance to TSS (bp)")

ggsave(file.path(opt$output_dir, "benchmark_median_distance_top50.png"), p2, width = 9, height = 6, dpi = 300)

p3 <- ggplot(summary_all[top_n == 50], aes(x = reorder(method, distal_frac), y = distal_frac, fill = method)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  theme_bw() +
  labs(title = sprintf("Distal fraction among top 50 links (> %s bp)", format(opt$distal_threshold, scientific = FALSE)),
       x = NULL, y = "Distal fraction")

ggsave(file.path(opt$output_dir, "benchmark_distal_fraction_top50.png"), p3, width = 9, height = 6, dpi = 300)

p4 <- ggplot(summary_all[top_n == 100], aes(x = reorder(method, unique_genes), y = unique_genes, fill = method)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  theme_bw() +
  labs(title = "Unique genes in top 100 links", x = NULL, y = "Unique genes")

ggsave(file.path(opt$output_dir, "benchmark_unique_genes_top100.png"), p4, width = 9, height = 6, dpi = 300)

# 2. Distance distributions of top links
plot_dt <- rbindlist(lapply(names(method_tables), function(nm) {
  x <- method_tables[[nm]][order(-score)][1:min(.N, opt$top_k_compare)]
  x[, method := nm]
  x
}), use.names = TRUE, fill = TRUE)
plot_dt <- plot_dt[is.finite(distance_bp)]

p5 <- ggplot(plot_dt, aes(x = distance_bp)) +
  geom_histogram(bins = 50) +
  facet_wrap(~method, scales = "free_y") +
  theme_bw() +
  labs(title = sprintf("Distance distribution of top %d links", opt$top_k_compare),
       x = "Distance to TSS (bp)", y = "Count")

ggsave(file.path(opt$output_dir, "benchmark_distance_distribution_topk.png"), p5, width = 14, height = 10, dpi = 300)

# 3. Overlap heatmap at top K pairs
ovk <- pair_overlap[k == opt$top_k_compare]
p6 <- ggplot(ovk, aes(x = method_a, y = method_b, fill = overlap_pairs)) +
  geom_tile() +
  geom_text(aes(label = overlap_pairs), size = 3) +
  theme_bw() +
  coord_equal() +
  labs(title = sprintf("Top-%d pair overlap across methods", opt$top_k_compare), x = NULL, y = NULL)

ggsave(file.path(opt$output_dir, "benchmark_pair_overlap_heatmap.png"), p6, width = 8, height = 7, dpi = 300)

# 4. LinkPeaks vs Reranker scatter on shared pairs
lr_shared <- merge(
  method_tables[["LinkPeaks"]][, .(peak, gene, linkpeaks = score)],
  method_tables[["Reranker"]][, .(peak, gene, reranker = score, distance_bp)],
  by = c("peak", "gene")
)
if (nrow(lr_shared) > 0) {
  p7 <- ggplot(lr_shared, aes(x = linkpeaks, y = reranker, color = distance_bp)) +
    geom_point(alpha = 0.5, size = 1) +
    theme_bw() +
    labs(title = "Shared-pair score comparison: LinkPeaks vs Reranker", x = "LinkPeaks score", y = "Reranker score")
  ggsave(file.path(opt$output_dir, "benchmark_linkpeaks_vs_reranker_scatter.png"), p7, width = 8, height = 6, dpi = 300)
}

# 5. Rank profile curves
rank_profile <- rbindlist(lapply(names(method_tables), function(nm) {
  x <- copy(method_tables[[nm]][order(-score)][1:min(.N, opt$top_k_compare)])
  x[, rank := seq_len(.N)]
  x[, method := nm]
  x
}), use.names = TRUE, fill = TRUE)

p8 <- ggplot(rank_profile, aes(x = rank, y = score, color = method)) +
  geom_line() +
  theme_bw() +
  labs(title = sprintf("Rank-score profiles across methods (top %d)", opt$top_k_compare), x = "Rank", y = "Method score")

ggsave(file.path(opt$output_dir, "benchmark_rank_score_profiles.png"), p8, width = 9, height = 6, dpi = 300)

# 6. Distal retention among top K
rank_profile[, distal := distance_bp > opt$distal_threshold]
distal_profile <- rank_profile[, .(distal_cumfrac = cumsum(distal) / seq_len(.N)), by = .(method, rank)]
p9 <- ggplot(distal_profile, aes(x = rank, y = distal_cumfrac, color = method)) +
  geom_line() +
  theme_bw() +
  labs(title = sprintf("Cumulative distal fraction through top %d ranks", opt$top_k_compare), x = "Rank", y = "Cumulative distal fraction")

ggsave(file.path(opt$output_dir, "benchmark_cumulative_distal_fraction.png"), p9, width = 9, height = 6, dpi = 300)

# ============================================================
# Console summary
# ============================================================
msg("Benchmark complete.")
msg("Methods compared: %s", paste(names(method_tables), collapse = ", "))
msg("Main summary: %s", file.path(opt$output_dir, "benchmark_summary_metrics.csv"))
msg("ORA summary: %s", file.path(opt$output_dir, "benchmark_ora_summary.csv"))
msg("Pair overlap: %s", file.path(opt$output_dir, "benchmark_pair_overlap.csv"))
msg("Standardized link tables written to: %s", opt$output_dir)
