#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(Signac)
  library(GenomeInfoDb)
  library(biovizBase)
  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
  library(TFBSTools)
  library(JASPAR2022)
  library(motifmatchr)
  library(SummarizedExperiment)
})

# ============================================================
# CLI
# ============================================================
option_list <- list(
  make_option(c("--data-dir"), type = "character", default = "data",
              help = "Directory containing 10x multiome input files [default %default]"),
  make_option(c("--h5-file"), type = "character", default = "filtered_feature_bc_matrix.h5",
              help = "10x h5 matrix filename inside --data-dir [default %default]"),
  make_option(c("--frag-file"), type = "character", default = "atac_fragments.tsv.gz",
              help = "ATAC fragments filename inside --data-dir [default %default]"),
  make_option(c("--candidate-top-k"), type = "integer", default = 5000,
              help = "Number of top LinkPeaks links to rerank [default %default]"),
  make_option(c("--link-distance"), type = "integer", default = 500000,
              help = "Maximum distance passed to LinkPeaks [default %default]"),
  make_option(c("--distance-d0"), type = "double", default = 50000,
              help = "Scale parameter for long-tail distance prior [default %default]"),
  make_option(c("--cluster-resolution"), type = "double", default = 0.5,
              help = "Resolution for FindClusters [default %default]"),
  make_option(c("--pca-dims"), type = "integer", default = 30,
              help = "Number of PCA dimensions for RNA [default %default]"),
  make_option(c("--lsi-dims-start"), type = "integer", default = 2,
              help = "Starting LSI dimension for ATAC [default %default]"),
  make_option(c("--lsi-dims-end"), type = "integer", default = 30,
              help = "Ending LSI dimension for ATAC [default %default]"),
  make_option(c("--species"), type = "integer", default = 9606,
              help = "NCBI taxonomy ID for JASPAR motifs [default %default]"),
  make_option(c("--collection"), type = "character", default = "CORE",
              help = "JASPAR collection [default %default]"),
  make_option(c("--tf-expressed-frac"), type = "double", default = 0.10,
              help = "Keep motifs whose TF is expressed in at least this fraction of cells [default %default]"),
  make_option(c("--top-k-motif-names"), type = "integer", default = 3,
              help = "How many top motif names to keep per peak [default %default]"),
  make_option(c("--output-dir"), type = "character", default = file.path("results", "pbmc", "features"),
              help = "Directory for feature-table outputs [default %default]"),
  make_option(c("--run-name"), type = "character", default = "pbmc",
              help = "Prefix for output files [default %default]"),
  make_option(c("--save-checkpoints"), action = "store_true", default = FALSE,
              help = "Save Seurat object checkpoints [default %default]"),
  make_option(c("--done-file"), type = "character", default = NULL,
              help = "Optional sentinel file touched on success"),
  make_option(c("--seed"), type = "integer", default = 42,
              help = "Random seed [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
names(opt) <- gsub("-", "_", names(opt), fixed = TRUE)
set.seed(opt$seed)

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

feature_file <- file.path(opt$output_dir, sprintf("%s_link_features.csv", opt$run_name))
baseline_file <- file.path(opt$output_dir, sprintf("%s_baseline_links_full.csv", opt$run_name))
baseline_dist_file <- file.path(opt$output_dir, sprintf("%s_baseline_links_with_distance.csv", opt$run_name))
checkpoint_file <- file.path(opt$output_dir, sprintf("%s_checkpoint_after_clustering.rds", opt$run_name))
post_linkpeaks_file <- file.path(opt$output_dir, sprintf("%s_post_linkpeaks_object.rds", opt$run_name))

# ============================================================
# Helpers
# ============================================================
msg <- function(...) cat(sprintf(...), "\n")

require_file <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
}

safe_scale <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

rescale01 <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) {
    out <- rep(0, length(x))
    out[is.na(x)] <- 0
    return(out)
  }
  out <- (x - rng[1]) / diff(rng)
  out[is.na(out)] <- 0
  out
}

safe_get_assay_data <- function(object, assay, what = "data") {
  tryCatch({
    GetAssayData(object, assay = assay, layer = what)
  }, error = function(e1) {
    tryCatch({
      GetAssayData(object, assay = assay, slot = what)
    }, error = function(e2) {
      stop(sprintf(
        "Could not extract assay data for assay='%s', what='%s'. layer-error: %s | slot-error: %s",
        assay, what, conditionMessage(e1), conditionMessage(e2)
      ))
    })
  })
}

normalize_tf_name <- function(x) {
  x <- toupper(as.character(x))
  x <- gsub("\\s+", "", x)
  x <- gsub("::", ";", x, fixed = TRUE)
  x <- gsub("\\(", ";", x)
  x <- gsub("\\)", "", x)
  x <- gsub("/", ";", x, fixed = TRUE)
  x <- gsub("-", "", x, fixed = TRUE)
  x
}

split_tf_tokens <- function(x) {
  toks <- unlist(strsplit(normalize_tf_name(x), ";", fixed = TRUE))
  unique(toks[nzchar(toks)])
}

build_tf_alias_map <- function(gene_names) {
  gene_names <- unique(as.character(gene_names))
  gene_names <- gene_names[!is.na(gene_names) & nzchar(gene_names)]
  alias_map <- setNames(vector("list", length(gene_names)), gene_names)

  for (g in gene_names) {
    toks <- unique(c(normalize_tf_name(g), gsub("-", "", toupper(g), fixed = TRUE)))
    alias_map[[g]] <- toks[nzchar(toks)]
  }

  alias_map
}

match_motifs_to_expressed_tfs <- function(motif_names, expressed_tfs) {
  alias_map <- build_tf_alias_map(expressed_tfs)
  motif_to_gene <- vector("list", length(motif_names))
  names(motif_to_gene) <- motif_names
  keep <- logical(length(motif_names))

  for (i in seq_along(motif_names)) {
    mtoks <- split_tf_tokens(motif_names[i])
    hits <- names(alias_map)[vapply(alias_map, function(atoks) any(mtoks %in% atoks), logical(1))]
    motif_to_gene[[i]] <- hits
    keep[i] <- length(hits) > 0
  }

  list(keep = keep, motif_to_gene = motif_to_gene)
}

parse_peak_table <- function(peak_vec) {
  peak_vec <- as.character(peak_vec)
  peak_vec2 <- gsub(":", "-", peak_vec, fixed = TRUE)
  parts <- strsplit(peak_vec2, "-", fixed = TRUE)
  ok <- lengths(parts) == 3

  if (!all(ok)) {
    bad <- unique(peak_vec[!ok])
    stop("Could not parse peak coordinates for: ", paste(head(bad, 5), collapse = ", "))
  }

  data.table(
    peak = peak_vec,
    peak_chr = vapply(parts, `[`, character(1), 1),
    peak_start = as.numeric(vapply(parts, `[`, character(1), 2)),
    peak_end = as.numeric(vapply(parts, `[`, character(1), 3))
  )
}

parse_peak_granges <- function(peak_vec) {
  peak_df <- parse_peak_table(peak_vec)
  gr <- GRanges(
    seqnames = peak_df$peak_chr,
    ranges = IRanges(start = peak_df$peak_start, end = peak_df$peak_end)
  )
  names(gr) <- peak_df$peak
  gr
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

extract_motif_score_matrix <- function(motif_obj, n_peaks, n_motifs) {
  mat <- tryCatch(as.matrix(motifScores(motif_obj)), error = function(e) NULL)

  if (is.null(mat)) {
    mat <- tryCatch(as.matrix(SummarizedExperiment::assay(motif_obj)), error = function(e) NULL)
  }

  if (is.null(mat)) stop("Could not extract continuous motif score matrix from motifmatchr result.")
  if (nrow(mat) == n_motifs && ncol(mat) == n_peaks) mat <- t(mat)
  if (nrow(mat) != n_peaks) stop("Unexpected motif score matrix dimensions after orientation check.")
  mat
}

get_top_names_from_scores <- function(score_vec, motif_names, k = 3) {
  score_vec <- as.numeric(score_vec)
  score_vec[!is.finite(score_vec)] <- NA_real_
  ord <- order(score_vec, decreasing = TRUE, na.last = NA)
  ord <- ord[score_vec[ord] > 0]
  if (length(ord) == 0) return(NA_character_)
  paste(unique(motif_names[head(ord, k)]), collapse = ";")
}

# ============================================================
# Input paths
# ============================================================
data_dir <- opt$data_dir
h5_file <- file.path(data_dir, opt$h5_file)
frag_file <- file.path(data_dir, opt$frag_file)
frag_index_file <- paste0(frag_file, ".tbi")

require_file(h5_file)
require_file(frag_file)
require_file(frag_index_file)

# ============================================================
# Build multiome object
# ============================================================
msg("Reading 10x multiome input...")
data <- Read10X_h5(h5_file)
msg("Matrices found: %s", paste(names(data), collapse = ", "))

if (!("Gene Expression" %in% names(data))) stop("Gene Expression matrix not found in h5 file.")
if (!("Peaks" %in% names(data))) stop("Peaks matrix not found in h5 file.")

msg("Building annotations...")
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"

msg("Creating Seurat object...")
obj <- CreateSeuratObject(counts = data$`Gene Expression`, assay = "RNA")
obj[["ATAC"]] <- CreateChromatinAssay(
  counts = data$Peaks,
  sep = c(":", "-"),
  genome = "hg38",
  fragments = frag_file,
  annotation = annotations
)

msg("Running RNA preprocessing...")
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, verbose = FALSE)

msg("Running ATAC preprocessing...")
DefaultAssay(obj) <- "ATAC"
obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = "q0")
obj <- RunSVD(obj)

msg("Running WNN integration...")
obj <- FindMultiModalNeighbors(
  obj,
  reduction.list = list("pca", "lsi"),
  dims.list = list(seq_len(opt$pca_dims), opt$lsi_dims_start:opt$lsi_dims_end)
)
obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", verbose = FALSE)
obj <- FindClusters(obj, graph.name = "wsnn", resolution = opt$cluster_resolution, verbose = FALSE)

if (isTRUE(opt$save_checkpoints)) {
  saveRDS(obj, checkpoint_file)
  msg("Saved checkpoint: %s", checkpoint_file)
}

# ============================================================
# LinkPeaks candidate generation
# ============================================================
msg("Running RegionStats + LinkPeaks...")
DefaultAssay(obj) <- "ATAC"
obj <- RegionStats(object = obj, genome = BSgenome.Hsapiens.UCSC.hg38)
obj <- LinkPeaks(
  object = obj,
  peak.assay = "ATAC",
  expression.assay = "RNA",
  distance = opt$link_distance
)

if (isTRUE(opt$save_checkpoints)) {
  saveRDS(obj, post_linkpeaks_file)
  msg("Saved post-LinkPeaks object: %s", post_linkpeaks_file)
}

links_df <- as.data.frame(Links(obj))
msg("Number of peak-gene links: %d", nrow(links_df))

baseline_df_full <- as.data.table(links_df[, c("peak", "gene", "score")])
baseline_df_full <- baseline_df_full[!is.na(score)]
setorder(baseline_df_full, peak, gene, -score)
baseline_df_full <- baseline_df_full[, .SD[1], by = .(peak, gene)]
setorder(baseline_df_full, -score)
fwrite(baseline_df_full, baseline_file)
msg("Wrote baseline links: %s", baseline_file)

# ============================================================
# Matrix access and coactivity features
# ============================================================
DefaultAssay(obj) <- "RNA"
rna_mat <- safe_get_assay_data(obj, assay = "RNA", what = "data")

DefaultAssay(obj) <- "ATAC"
atac_mat <- safe_get_assay_data(obj, assay = "ATAC", what = "data")

msg("Scoring top %d baseline links...", opt$candidate_top_k)
candidates <- head(baseline_df_full, opt$candidate_top_k)

candidate_genes <- intersect(unique(candidates$gene), rownames(rna_mat))
candidate_peaks <- intersect(unique(candidates$peak), rownames(atac_mat))

rna_list <- setNames(lapply(candidate_genes, function(g) safe_scale(as.numeric(rna_mat[g, ]))), candidate_genes)
atac_list <- setNames(lapply(candidate_peaks, function(p) safe_scale(as.numeric(atac_mat[p, ]))), candidate_peaks)

score_rows <- vector("list", nrow(candidates))
keep_idx <- logical(nrow(candidates))

for (i in seq_len(nrow(candidates))) {
  peak <- candidates$peak[i]
  gene <- candidates$gene[i]
  if (!(gene %in% names(rna_list))) next
  if (!(peak %in% names(atac_list))) next

  rna_z <- rna_list[[gene]]
  atac_z <- atac_list[[peak]]

  rna_bin <- rna_z > 1
  atac_bin <- atac_z > 1
  activity_penalty <- mean(rna_bin, na.rm = TRUE) * mean(atac_bin, na.rm = TRUE)

  keep_idx[i] <- TRUE
  score_rows[[i]] <- data.frame(
    peak = peak,
    gene = gene,
    link_score = candidates$score[i],
    add = mean(pmax(rna_z, 0) + pmax(atac_z, 0), na.rm = TRUE),
    mul = mean(rna_z * atac_z, na.rm = TRUE),
    mul_weigh = mean(pmax(rna_z, 0) * pmax(atac_z, 0), na.rm = TRUE),
    mul_strict = mean(rna_bin & atac_bin, na.rm = TRUE),
    adj = mean(rna_bin & atac_bin, na.rm = TRUE) / (activity_penalty + 1e-6),
    stringsAsFactors = FALSE
  )
}

features <- rbindlist(score_rows[keep_idx], fill = TRUE)
if (nrow(features) == 0) stop("No candidate links survived scoring. Check feature names and matrices.")
msg("Scored candidate links retained: %d", nrow(features))

# ============================================================
# Distance features
# ============================================================
msg("Adding transcript-derived TSS distance features...")
gene_coord_df <- build_gene_tss_table()

peak_df <- unique(parse_peak_table(features$peak))
peak_df[, peak_mid := (peak_start + peak_end) / 2]

features <- merge(features, peak_df, by = "peak", all.x = TRUE, sort = FALSE)
features <- merge(features, gene_coord_df[, .(gene, gene_chr, tss, tx_id, tx_biotype)], by = "gene", all.x = TRUE, sort = FALSE)
features[, distance_bp := ifelse(is.na(gene_chr) | peak_chr != gene_chr, Inf, abs(peak_mid - tss))]
features[, distance_score := ifelse(is.finite(distance_bp), 1 / (1 + (distance_bp / opt$distance_d0)^2), 0)]

baseline_peak_df <- unique(parse_peak_table(baseline_df_full$peak))
baseline_peak_df[, peak_mid := (peak_start + peak_end) / 2]
baseline_dist <- merge(baseline_df_full, baseline_peak_df, by = "peak", all.x = TRUE, sort = FALSE)
baseline_dist <- merge(baseline_dist, gene_coord_df[, .(gene, gene_chr, tss)], by = "gene", all.x = TRUE, sort = FALSE)
baseline_dist[, distance_bp := ifelse(is.na(gene_chr) | peak_chr != gene_chr, Inf, abs(peak_mid - tss))]
baseline_dist[, distance_score := ifelse(is.finite(distance_bp), 1 / (1 + (distance_bp / opt$distance_d0)^2), 0)]
fwrite(baseline_dist, baseline_dist_file)
msg("Wrote baseline distance table: %s", baseline_dist_file)

# ============================================================
# Continuous motif + TF features
# ============================================================
msg("Running motif scoring with continuous motif support...")
peak_ids <- unique(features$peak)
peak_gr <- parse_peak_granges(peak_ids)
common_seqlevels <- intersect(seqlevels(peak_gr), seqlevels(BSgenome.Hsapiens.UCSC.hg38))
peak_gr <- keepSeqlevels(peak_gr, common_seqlevels, pruning.mode = "coarse")
peak_gr <- peak_gr[!is.na(start(peak_gr)) & !is.na(end(peak_gr)) & width(peak_gr) > 0]
msg("Peaks retained for motif scoring: %d", length(peak_gr))

jaspar_opts <- list(collection = opt$collection, tax_group = "vertebrates", species = opt$species)
pfm_list <- getMatrixSet(JASPAR2022, jaspar_opts)
if (length(pfm_list) == 0) stop("No JASPAR motifs returned. Check species/collection settings.")

motif_names <- vapply(pfm_list, name, character(1))
expressed_tfs <- rownames(rna_mat)[Matrix::rowMeans(rna_mat > 0) > opt$tf_expressed_frac]
tf_match <- match_motifs_to_expressed_tfs(motif_names = motif_names, expressed_tfs = expressed_tfs)

pfm_list <- pfm_list[tf_match$keep]
motif_names <- motif_names[tf_match$keep]
motif_to_gene <- tf_match$motif_to_gene[tf_match$keep]
if (length(pfm_list) == 0) stop("No JASPAR motifs left after TF-expression matching.")
msg("Motifs retained after TF expression filtering: %d", length(pfm_list))

motif_match_obj <- matchMotifs(
  pwms = pfm_list,
  subject = peak_gr,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  out = "scores",
  bg = "genome"
)

motif_score_mat <- extract_motif_score_matrix(motif_match_obj, n_peaks = length(peak_gr), n_motifs = length(pfm_list))
rownames(motif_score_mat) <- names(peak_gr)
colnames(motif_score_mat) <- motif_names

if (ncol(motif_score_mat) > 0) {
  scaled_cols <- lapply(seq_len(ncol(motif_score_mat)), function(j) rescale01(motif_score_mat[, j]))
  motif_score_mat <- matrix(
    do.call(cbind, scaled_cols),
    nrow = length(peak_gr),
    ncol = length(motif_names),
    dimnames = list(names(peak_gr), motif_names)
  )
}

gene_mean_expr <- Matrix::rowMeans(rna_mat)
motif_tf_expr <- vapply(seq_along(motif_to_gene), function(i) {
  genes <- intersect(motif_to_gene[[i]], names(gene_mean_expr))
  if (length(genes) == 0) return(0)
  max(gene_mean_expr[genes], na.rm = TRUE)
}, numeric(1))

names(motif_tf_expr) <- colnames(motif_score_mat)
motif_tf_expr <- rescale01(motif_tf_expr)

peak_tf_score <- as.numeric(motif_score_mat %*% motif_tf_expr)
peak_tf_score <- rescale01(peak_tf_score)
names(peak_tf_score) <- rownames(motif_score_mat)

peak_best_motif_score <- apply(motif_score_mat, 1, max, na.rm = TRUE)
peak_best_motif_score[!is.finite(peak_best_motif_score)] <- 0
peak_best_motif_score <- rescale01(peak_best_motif_score)

peak_top_names <- apply(
  motif_score_mat,
  1,
  function(x) get_top_names_from_scores(x, colnames(motif_score_mat), k = opt$top_k_motif_names)
)

peak_summary <- data.table(
  peak = rownames(motif_score_mat),
  motif_score = unname(peak_best_motif_score),
  tf_score = unname(peak_tf_score),
  motif_names = unname(peak_top_names)
)

features <- merge(features, peak_summary, by = "peak", all.x = TRUE, sort = FALSE)
features[is.na(motif_score), motif_score := 0]
features[is.na(tf_score), tf_score := 0]
features[is.na(motif_names), motif_names := NA_character_]

setorder(features, peak, gene, -link_score)
features <- features[, .SD[1], by = .(peak, gene)]
setorder(features, -link_score)

fwrite(features, feature_file)
msg("Wrote feature table: %s", feature_file)

if (!is.null(opt$done_file)) {
  dir.create(dirname(opt$done_file), showWarnings = FALSE, recursive = TRUE)
  file.create(opt$done_file)
  msg("Touched done file: %s", opt$done_file)
}

msg("Done.")
