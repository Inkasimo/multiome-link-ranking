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
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
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
  make_option(c("--lambda-distance"), type = "double", default = 0.30,
              help = "Strength of distance modifier in final score [default %default]"),
  make_option(c("--alpha-tf"), type = "double", default = 0.30,
              help = "Strength of TF modifier in final score [default %default]"),
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
  make_option(c("--motif-min-score"), type = "character", default = NULL,
            help = "Unused in continuous-score mode; kept for backward compatibility [default %default]"),
  make_option(c("--tf-expressed-frac"), type = "double", default = 0.10,
              help = "Keep motifs whose TF is expressed in at least this fraction of cells [default %default]"),
  make_option(c("--top-k-motif-names"), type = "integer", default = 3,
              help = "How many top motif names to keep per peak [default %default]"),
  make_option(c("--ora-top-n"), type = "integer", default = 100,
              help = "Top N genes per method for ORA [default %default]"),
  make_option(c("--ora-show-category"), type = "integer", default = 15,
              help = "Number of categories to show in ORA dotplots [default %default]"),
  make_option(c("--tier-high-quantile"), type = "double", default = 0.90,
              help = "Quantile cutoff for High tier [default %default]"),
  make_option(c("--tier-medium-quantile"), type = "double", default = 0.70,
              help = "Quantile cutoff for Medium tier [default %default]"),
  make_option(c("--output-prefix"), type = "character", default = "multiome_rie",
              help = "Prefix for output files [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ============================================================
# Helpers
# ============================================================
msg <- function(...) cat(sprintf(...), "\n")

require_file <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }
}

safe_scale <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  }
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
  out <- tryCatch({
    GetAssayData(object, assay = assay, layer = what)
  }, error = function(e1) {
    tryCatch({
      GetAssayData(object, assay = assay, slot = what)
    }, error = function(e2) {
      stop(
        sprintf(
          "Could not extract assay data for assay='%s', what='%s'. layer-error: %s | slot-error: %s",
          assay, what, conditionMessage(e1), conditionMessage(e2)
        )
      )
    })
  })
  out
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
  x <- normalize_tf_name(x)
  toks <- unlist(strsplit(x, ";", fixed = TRUE))
  toks <- toks[nzchar(toks)]
  unique(toks)
}

build_tf_alias_map <- function(gene_names) {
  gene_names <- unique(as.character(gene_names))
  gene_names <- gene_names[!is.na(gene_names) & nzchar(gene_names)]

  alias_map <- setNames(vector("list", length(gene_names)), gene_names)

  for (g in gene_names) {
    toks <- unique(c(
      normalize_tf_name(g),
      gsub("-", "", toupper(g), fixed = TRUE)
    ))
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

    hits <- names(alias_map)[vapply(alias_map, function(atoks) {
      any(mtoks %in% atoks)
    }, logical(1))]

    motif_to_gene[[i]] <- hits
    keep[i] <- length(hits) > 0
  }

  list(
    keep = keep,
    motif_to_gene = motif_to_gene
  )
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

  data.frame(
    peak = peak_vec,
    peak_chr = vapply(parts, `[`, character(1), 1),
    peak_start = as.numeric(vapply(parts, `[`, character(1), 2)),
    peak_end = as.numeric(vapply(parts, `[`, character(1), 3)),
    stringsAsFactors = FALSE
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

build_gene_coord_table <- function(annotations) {
  gene_names <- mcols(annotations)$gene_name
  seqs <- as.character(seqnames(annotations))
  starts <- start(annotations)
  ends <- end(annotations)
  strands <- as.character(strand(annotations))

  keep <- !is.na(gene_names) & nzchar(gene_names)

  ann_df <- data.frame(
    gene = gene_names[keep],
    gene_chr = seqs[keep],
    gene_start = starts[keep],
    gene_end = ends[keep],
    gene_strand = strands[keep],
    stringsAsFactors = FALSE
  )

  ann_df$tss <- ifelse(ann_df$gene_strand == "-", ann_df$gene_end, ann_df$gene_start)

  gene_coord_df <- aggregate(
    tss ~ gene + gene_chr,
    data = ann_df,
    FUN = function(x) {
      ux <- unique(x)
      ux[which.min(abs(ux - stats::median(ux)))]
    }
  )

  gene_coord_df
}

extract_motif_score_matrix <- function(motif_obj, n_peaks, n_motifs) {
  # Try common extractors across motifmatchr / SummarizedExperiment versions.
  mat <- NULL

  mat <- tryCatch({
    as.matrix(motifScores(motif_obj))
  }, error = function(e) NULL)

  if (is.null(mat)) {
    mat <- tryCatch({
      as.matrix(SummarizedExperiment::assay(motif_obj))
    }, error = function(e) NULL)
  }

  if (is.null(mat)) {
    stop("Could not extract continuous motif score matrix from motifmatchr result.")
  }

  if (nrow(mat) == n_motifs && ncol(mat) == n_peaks) {
    mat <- t(mat)
  }

  if (nrow(mat) != n_peaks) {
    stop("Unexpected motif score matrix dimensions after orientation check.")
  }

  mat
}

get_top_names_from_scores <- function(score_vec, motif_names, k = 3) {
  score_vec <- as.numeric(score_vec)
  score_vec[!is.finite(score_vec)] <- NA_real_
  ord <- order(score_vec, decreasing = TRUE, na.last = NA)
  ord <- ord[score_vec[ord] > 0]
  if (length(ord) == 0) {
    return(NA_character_)
  }
  paste(unique(motif_names[head(ord, k)]), collapse = ";")
}

assign_tiers <- function(score_vec, high_q = 0.90, medium_q = 0.70) {
  if (medium_q >= high_q) {
    stop("tier-medium-quantile must be < tier-high-quantile")
  }
  hi <- as.numeric(stats::quantile(score_vec, probs = high_q, na.rm = TRUE))
  med <- as.numeric(stats::quantile(score_vec, probs = medium_q, na.rm = TRUE))
  out <- ifelse(score_vec >= hi, "High",
                ifelse(score_vec >= med, "Medium", "Low"))
  factor(out, levels = c("High", "Medium", "Low"))
}

safe_bitr <- function(genes) {
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (length(genes) == 0) {
    return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
  }
  out <- suppressWarnings(
    bitr(
      genes,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
  )
  unique(out)
}

safe_enrich_go <- function(entrez_ids, universe_ids) {
  if (length(entrez_ids) == 0 || length(universe_ids) == 0) {
    return(NULL)
  }
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
  }, error = function(e) {
    warning("enrichGO failed: ", conditionMessage(e))
    NULL
  })
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

# ============================================================
# Input paths
# ============================================================
data_dir <- opt$data_dir
h5_file <- file.path(data_dir, opt$h5_file)
frag_file <- file.path(data_dir, opt$frag_file)

require_file(h5_file)
require_file(frag_file)

msg("Reading 10x multiome input...")
data <- Read10X_h5(h5_file)
msg("Matrices found: %s", paste(names(data), collapse = ", "))

if (!("Gene Expression" %in% names(data))) {
  stop("Gene Expression matrix not found in h5 file.")
}
if (!("Peaks" %in% names(data))) {
  stop("Peaks matrix not found in h5 file.")
}

# ============================================================
# Build multiome object
# ============================================================
msg("Building annotations...")
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"

automatic_prefix <- opt$output_prefix
obj_checkpoint <- sprintf("%s_checkpoint_after_clustering.rds", automatic_prefix)
obj_final_rds <- sprintf("%s_multiome_object.rds", automatic_prefix)

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

saveRDS(obj, obj_checkpoint)
msg("Saved checkpoint: %s", obj_checkpoint)

# ============================================================
# Baseline candidate generation with LinkPeaks
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

links <- Links(obj)
links_df <- as.data.frame(links)
msg("Number of peak-gene links: %d", nrow(links_df))

saveRDS(obj, obj_final_rds)
msg("Saved full multiome object: %s", obj_final_rds)

baseline_df <- links_df[, c("peak", "gene", "score")]
baseline_df <- baseline_df[!is.na(baseline_df$score), ]
baseline_df <- baseline_df[order(-baseline_df$score), ]

baseline_df <- as.data.table(baseline_df)
setorder(baseline_df, peak, gene, -score)
baseline_df <- baseline_df[, .SD[1], by = .(peak, gene)]
setorder(baseline_df, -score)

write.csv(baseline_df, sprintf("%s_baseline_links.csv", automatic_prefix), row.names = FALSE)

# ============================================================
# Access matrices
# ============================================================
DefaultAssay(obj) <- "RNA"
rna_mat <- safe_get_assay_data(obj, assay = "RNA", what = "data")

DefaultAssay(obj) <- "ATAC"
atac_mat <- safe_get_assay_data(obj, assay = "ATAC", what = "data")

# ============================================================
# Core reranking: additive, multiplicative, strict, adjusted
# ============================================================
msg("Scoring top %d baseline links...", opt$candidate_top_k)
candidates <- head(baseline_df, opt$candidate_top_k)

score_rows <- vector("list", nrow(candidates))
keep_idx <- logical(nrow(candidates))

for (i in seq_len(nrow(candidates))) {
  peak <- candidates$peak[i]
  gene <- candidates$gene[i]

  if (!(gene %in% rownames(rna_mat))) next
  if (!(peak %in% rownames(atac_mat))) next

  rna <- as.numeric(rna_mat[gene, ])
  atac <- as.numeric(atac_mat[peak, ])

  rna_z <- safe_scale(rna)
  atac_z <- safe_scale(atac)

  link_score <- candidates$score[i]
  score_add <- mean(pmax(rna_z, 0) + pmax(atac_z, 0), na.rm = TRUE)
  score_mul <- mean(rna_z * atac_z, na.rm = TRUE)
  score_mul_weighted <- mean(pmax(rna_z, 0) * pmax(atac_z, 0), na.rm = TRUE)

  rna_bin <- rna_z > 1
  atac_bin <- atac_z > 1
  score_mul_strict <- mean(rna_bin & atac_bin, na.rm = TRUE)
  activity_penalty <- mean(rna_bin, na.rm = TRUE) * mean(atac_bin, na.rm = TRUE)
  score_adj <- score_mul_strict / (activity_penalty + 1e-6)

  keep_idx[i] <- TRUE
  score_rows[[i]] <- data.frame(
    peak = peak,
    gene = gene,
    link_score = link_score,
    add = score_add,
    mul = score_mul,
    mul_weigh = score_mul_weighted,
    mul_strict = score_mul_strict,
    adj = score_adj,
    stringsAsFactors = FALSE
  )
}

results <- data.table::rbindlist(score_rows[keep_idx], fill = TRUE)
if (nrow(results) == 0) {
  stop("No candidate links survived scoring. Check feature names and matrices.")
}

write.csv(results, sprintf("%s_test_scores.csv", automatic_prefix), row.names = FALSE)
msg("Scored candidate links retained: %d", nrow(results))

# ============================================================
# Distance prior
# ============================================================
msg("Adding distance prior...")
gene_coord_df <- build_gene_coord_table(annotations)
peak_df <- parse_peak_table(results$peak)
peak_df$peak_mid <- (peak_df$peak_start + peak_df$peak_end) / 2

results <- merge(results, peak_df, by = "peak", all.x = TRUE, sort = FALSE)
results <- merge(
  results,
  gene_coord_df[, c("gene", "gene_chr", "tss")],
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)

results[, distance_bp := ifelse(
  is.na(gene_chr) | peak_chr != gene_chr,
  Inf,
  abs(peak_mid - tss)
)]

results[, distance_score := ifelse(
  is.finite(distance_bp),
  1 / (1 + (distance_bp / opt$distance_d0)^2),
  0
)]

results[, final_v5 := mul_weigh * ((1 - opt$lambda_distance) + opt$lambda_distance * distance_score)]
results[, rank_link := frank(-link_score, ties.method = "average")]
results[, rank_final_v5 := frank(-final_v5, ties.method = "average")]
results[, rank_diff_v5 := rank_link - rank_final_v5]

# ============================================================
# Continuous motif + TF layer
# ============================================================
msg("Running motif scoring with continuous motif support...")
peak_ids <- unique(results$peak)
peak_gr <- parse_peak_granges(peak_ids)
common_seqlevels <- intersect(seqlevels(peak_gr), seqlevels(BSgenome.Hsapiens.UCSC.hg38))
peak_gr <- keepSeqlevels(peak_gr, common_seqlevels, pruning.mode = "coarse")
peak_gr <- peak_gr[!is.na(start(peak_gr)) & !is.na(end(peak_gr)) & width(peak_gr) > 0]

msg("Peaks retained for motif scoring: %d", length(peak_gr))

jaspar_opts <- list(
  collection = opt$collection,
  tax_group = "vertebrates",
  species = opt$species
)

pfm_list <- getMatrixSet(JASPAR2022, jaspar_opts)
if (length(pfm_list) == 0) {
  stop("No JASPAR motifs returned. Check species/collection settings.")
}

motif_names <- vapply(pfm_list, name, character(1))
expressed_tfs <- rownames(rna_mat)[Matrix::rowMeans(rna_mat > 0) > opt$tf_expressed_frac]

tf_match <- match_motifs_to_expressed_tfs(
  motif_names = motif_names,
  expressed_tfs = expressed_tfs
)

pfm_list <- pfm_list[tf_match$keep]
motif_names <- motif_names[tf_match$keep]
motif_to_gene <- tf_match$motif_to_gene[tf_match$keep]

if (length(pfm_list) == 0) {
  stop("No JASPAR motifs left after TF-expression matching.")
}

msg("Motifs retained after TF expression filtering: %d", length(pfm_list))

motif_match_obj <- matchMotifs(
  pwms = pfm_list,
  subject = peak_gr,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  out = "scores",
  bg = "genome"
)

motif_score_mat <- extract_motif_score_matrix(
  motif_obj = motif_match_obj,
  n_peaks = length(peak_gr),
  n_motifs = length(pfm_list)
)
rownames(motif_score_mat) <- names(peak_gr)
colnames(motif_score_mat) <- motif_names

# Rescale each motif column to 0-1 to reduce motif-specific scale effects.
if (ncol(motif_score_mat) > 0) {
  motif_score_mat <- apply(motif_score_mat, 2, rescale01)
  motif_score_mat <- as.matrix(motif_score_mat)
  if (nrow(motif_score_mat) != length(peak_gr)) {
    motif_score_mat <- t(motif_score_mat)
  }
  rownames(motif_score_mat) <- names(peak_gr)
  colnames(motif_score_mat) <- motif_names
}

# TF expression vector aligned to motif columns.
gene_mean_expr <- Matrix::rowMeans(rna_mat)

motif_tf_expr <- vapply(seq_along(motif_to_gene), function(i) {
  genes <- motif_to_gene[[i]]
  genes <- intersect(genes, names(gene_mean_expr))
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
  motif_tf_score = unname(peak_tf_score),
  motif_names = unname(peak_top_names)
)

results <- merge(results, peak_summary, by = "peak", all.x = TRUE, sort = FALSE)
results[is.na(motif_score), motif_score := 0]
results[is.na(tf_score), tf_score := 0]
results[is.na(motif_tf_score), motif_tf_score := 0]

# compute final_v6 first
results[, final_v6 := final_v5 * (1 + opt$alpha_tf * tf_score)]

# deduplicate on peak-gene
results <- as.data.table(results)
setorder(results, peak, gene, -final_v6, -link_score)
results <- results[, .SD[1], by = .(peak, gene)]

# recompute ranks after dedup
results[, rank_link := frank(-link_score, ties.method = "average")]
results[, rank_final_v5 := frank(-final_v5, ties.method = "average")]
results[, rank_diff_v5 := rank_link - rank_final_v5]
results[, rank_final_v6 := frank(-final_v6, ties.method = "average")]
results[, rank_diff_v6 := rank_link - rank_final_v6]
results[, tier := assign_tiers(final_v6, opt$tier_high_quantile, opt$tier_medium_quantile)]
setorder(results, -final_v6)
write.csv(results, sprintf("%s_ranked_links.csv", automatic_prefix), row.names = FALSE)

# Additional summary tables
write.csv(
  head(results[, c("gene", "peak", "link_score", "mul_weigh", "distance_bp", "distance_score", "tf_score", "final_v6", "tier")], 100),
  sprintf("%s_top100_final_links.csv", automatic_prefix),
  row.names = FALSE
)

promoted <- results[order(-rank_diff_v6), c("gene", "peak", "link_score", "final_v6", "distance_bp", "distance_score", "tf_score", "rank_link", "rank_final_v6", "rank_diff_v6")]
demoted <- results[order(rank_diff_v6), c("gene", "peak", "link_score", "final_v6", "distance_bp", "distance_score", "tf_score", "rank_link", "rank_final_v6", "rank_diff_v6")]
write.csv(head(promoted, 100), sprintf("%s_top_promoted_links.csv", automatic_prefix), row.names = FALSE)
write.csv(head(demoted, 100), sprintf("%s_top_demoted_links.csv", automatic_prefix), row.names = FALSE)

# ============================================================
# Tiers summary
# ============================================================
tier_summary <- results[, .(
  n_links = .N,
  unique_genes = uniqueN(gene),
  median_final_v6 = median(final_v6, na.rm = TRUE),
  median_distance_bp = median(distance_bp[is.finite(distance_bp)], na.rm = TRUE)
), by = tier]
write.csv(tier_summary, sprintf("%s_tier_summary.csv", automatic_prefix), row.names = FALSE)

# ============================================================
# ORA: baseline vs final combined ranking
# ============================================================
msg("Running ORA on top %d genes...", opt$ora_top_n)
baseline_genes <- unique(head(results[order(-link_score)]$gene, opt$ora_top_n))
final_genes <- unique(head(results[order(-final_v6)]$gene, opt$ora_top_n))
background_genes <- unique(results$gene)

baseline_map <- safe_bitr(baseline_genes)
final_map <- safe_bitr(final_genes)
background_map <- safe_bitr(background_genes)

baseline_entrez <- unique(baseline_map$ENTREZID)
final_entrez <- unique(final_map$ENTREZID)
background_entrez <- unique(background_map$ENTREZID)

ora_baseline <- safe_enrich_go(baseline_entrez, background_entrez)
ora_final <- safe_enrich_go(final_entrez, background_entrez)

baseline_ora_df <- if (is.null(ora_baseline)) data.frame() else as.data.frame(ora_baseline)
final_ora_df <- if (is.null(ora_final)) data.frame() else as.data.frame(ora_final)

write.csv(baseline_ora_df, sprintf("%s_ora_baseline_GO_BP.csv", automatic_prefix), row.names = FALSE)
write.csv(final_ora_df, sprintf("%s_ora_final_v6_GO_BP.csv", automatic_prefix), row.names = FALSE)

baseline_dot <- make_dotplot_safe(ora_baseline, "ORA: Baseline LinkPeaks (top genes)", opt$ora_show_category)
final_dot <- make_dotplot_safe(ora_final, "ORA: final_v6 (top genes)", opt$ora_show_category)

ggsave(sprintf("%s_ora_baseline_dotplot.png", automatic_prefix), baseline_dot, width = 10, height = 7, dpi = 300)
ggsave(sprintf("%s_ora_final_v6_dotplot.png", automatic_prefix), final_dot, width = 10, height = 7, dpi = 300)

if (nrow(baseline_ora_df) > 0 || nrow(final_ora_df) > 0) {
  baseline_top <- head(baseline_ora_df, 10)
  final_top <- head(final_ora_df, 10)
  if (nrow(baseline_top) > 0) baseline_top$method <- "baseline"
  if (nrow(final_top) > 0) final_top$method <- "final_v6"

  ora_compare <- rbind(
    if (nrow(baseline_top) > 0) baseline_top[, c("Description", "Count", "p.adjust", "method")] else NULL,
    if (nrow(final_top) > 0) final_top[, c("Description", "Count", "p.adjust", "method")] else NULL
  )

  if (!is.null(ora_compare) && nrow(ora_compare) > 0) {
    ora_compare$log10_padj <- -log10(ora_compare$p.adjust)
    ora_compare$Description <- factor(ora_compare$Description, levels = rev(unique(ora_compare$Description)))

    p_compare <- ggplot(ora_compare, aes(x = Count, y = Description, color = method, size = log10_padj)) +
      geom_point() +
      facet_wrap(~method, scales = "free_y") +
      theme_bw() +
      ggtitle("GO BP ORA: baseline vs final_v6")
  } else {
    p_compare <- empty_plot("GO BP ORA: baseline vs final_v6")
  }
} else {
  p_compare <- empty_plot("GO BP ORA: baseline vs final_v6")
}

ggsave(sprintf("%s_ora_baseline_vs_final_v6.png", automatic_prefix), p_compare, width = 12, height = 8, dpi = 300)

# ============================================================
# Ranking summary metrics
# ============================================================
summary_dt <- data.table(
  metric = c(
    "n_candidate_links_scored",
    "cor_link_vs_mul_weigh",
    "cor_link_vs_final_v5",
    "cor_link_vs_final_v6",
    "median_distance_top50_linkpeaks",
    "median_distance_top50_final_v6",
    "distal_frac_top50_final_v6_gt50kb",
    "n_unique_baseline_genes_topN",
    "n_unique_final_genes_topN",
    "n_baseline_ora_terms",
    "n_final_ora_terms"
  ),
  value = c(
    nrow(results),
    suppressWarnings(cor(results$link_score, results$mul_weigh, use = "complete.obs")),
    suppressWarnings(cor(results$link_score, results$final_v5, use = "complete.obs")),
    suppressWarnings(cor(results$link_score, results$final_v6, use = "complete.obs")),
    median(head(results[order(-link_score)]$distance_bp, 50), na.rm = TRUE),
    median(head(results[order(-final_v6)]$distance_bp, 50), na.rm = TRUE),
    mean(head(results[order(-final_v6)]$distance_bp, 50) > 50000, na.rm = TRUE),
    length(unique(baseline_genes)),
    length(unique(final_genes)),
    nrow(baseline_ora_df),
    nrow(final_ora_df)
  )
)
write.csv(summary_dt, sprintf("%s_summary_metrics.csv", automatic_prefix), row.names = FALSE)

# Top-N overlap summary
noverlap <- lapply(c(10, 20, 50, 100, 200), function(n) {
  top_link_n <- unique(head(results[order(-link_score), .(peak, gene)], n))
  top_final_n <- unique(head(results[order(-final_v6), .(peak, gene)], n))

  overlap_n <- nrow(merge(
    as.data.frame(top_link_n),
    as.data.frame(top_final_n),
    by = c("peak", "gene")
  ))

  data.frame(top_n = n, overlap = overlap_n)
})
noverlap <- do.call(rbind, noverlap)
write.csv(noverlap, sprintf("%s_topN_overlap_baseline_vs_final_v6.csv", automatic_prefix), row.names = FALSE)

# Scatter plot baseline vs final score
p_scatter <- ggplot(results, aes(x = link_score, y = final_v6, color = tier)) +
  geom_point(alpha = 0.5, size = 1) +
  theme_bw() +
  ggtitle("Baseline LinkPeaks vs final_v6") +
  xlab("Baseline LinkPeaks score") +
  ylab("final_v6 score")
ggsave(sprintf("%s_baseline_vs_final_v6_scatter.png", automatic_prefix), p_scatter, width = 8, height = 6, dpi = 300)

# Distance distribution among top-ranked links
results$top200_final_v6 <- seq_len(nrow(results)) <= 200
p_dist <- ggplot(results[is.finite(distance_bp), ], aes(x = distance_bp, fill = top200_final_v6)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  theme_bw() +
  ggtitle("Distance distribution: all links vs top200 final_v6") +
  xlab("Distance to gene TSS (bp)") +
  ylab("Count")
ggsave(sprintf("%s_distance_distribution.png", automatic_prefix), p_dist, width = 8, height = 6, dpi = 300)

# ============================================================
# Console summary
# ============================================================
msg("Done.")
msg("Top baseline genes (n=%d unique): %d", opt$ora_top_n, length(unique(baseline_genes)))
msg("Top final_v6 genes (n=%d unique): %d", opt$ora_top_n, length(unique(final_genes)))
msg("Tier counts: %s", paste(capture.output(print(table(results$tier))), collapse = " "))
msg("Main output: %s_ranked_links.csv", automatic_prefix)
