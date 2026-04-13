library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(biovizBase)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

data_dir <- "data/"
h5_file  <- file.path(data_dir, "filtered_feature_bc_matrix.h5")
frag_file <- file.path(data_dir, "atac_fragments.tsv.gz")

data <- Read10X_h5(h5_file)
print(names(data))

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"

obj <- CreateSeuratObject(
  counts = data$`Gene Expression`,
  assay = "RNA"
)

obj[["ATAC"]] <- CreateChromatinAssay(
  counts = data$Peaks,
  sep = c(":", "-"),
  genome = "hg38",
  fragments = frag_file,
  annotation = annotations
)


DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)

DefaultAssay(obj) <- "ATAC"
obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj)
obj <- RunSVD(obj)

obj <- FindMultiModalNeighbors(
  obj,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:30, 2:30)
)

obj <- RunUMAP(
  obj,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap"
)

obj <- FindClusters(
  obj,
  graph.name = "wsnn",
  resolution = 0.5
)

print(DimPlot(obj, reduction = "wnn.umap", label = TRUE))
print(DimPlot(obj, reduction = "pca"))
print(DimPlot(obj, reduction = "lsi"))

saveRDS(obj, "checkpoint_after_clustering.rds")


DefaultAssay(obj) <- "ATAC"

obj <- RegionStats(
  object = obj,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

obj <- LinkPeaks(
  object = obj,
  peak.assay = "ATAC",
  expression.assay = "RNA",
  distance = 500000
)

links <- Links(obj)
cat("Number of peak-gene links:", nrow(as.data.frame(links)), "\n")

saveRDS(obj, "pilot_multiome.rds")

##########################
df <- as.data.frame(links)

# keep baseline output
df <- df[, c("peak", "gene", "score")]
df <- df[!is.na(df$score), ]
df <- df[order(-df$score), ]

write.csv(df, "baseline_links.csv", row.names = FALSE)

# expression and accessibility matrices
DefaultAssay(obj) <- "RNA"
rna_mat <- GetAssayData(obj, layer = "data")

DefaultAssay(obj) <- "ATAC"
atac_mat <- GetAssayData(obj, layer = "data")

# use top baseline links as pilot candidates
df_small <- head(df, 2000)

results <- data.frame()

for (i in 1:nrow(df_small)) {
  peak <- df_small$peak[i]
  gene <- df_small$gene[i]

  if (!(gene %in% rownames(rna_mat))) next
  if (!(peak %in% rownames(atac_mat))) next

  rna <- as.numeric(rna_mat[gene, ])
  atac <- as.numeric(atac_mat[peak, ])

  # normalize per link so one modality doesn't dominate by scale alone
  rna_z <- as.numeric(scale(rna))
  atac_z <- as.numeric(scale(atac))

  # baseline from LinkPeaks
  link_score <- df_small$score[i]

  # alternative scoring rules
  score_add <- mean(pmax(rna_z, 0) + pmax(atac_z, 0), na.rm = TRUE)
  score_mul <- mean(rna_z * atac_z, na.rm = TRUE)
  score_mul_weighted <- mean(pmax(rna_z, 0) * pmax(atac_z, 0), na.rm = TRUE)

  # stricter co-activation
  rna_bin <- rna_z > 1
  atac_bin <- atac_z > 1
  score_mul_strict <- mean(rna_bin & atac_bin, na.rm = TRUE)

  # penalize globally active links
  activity_penalty <- mean(rna_bin, na.rm = TRUE) * mean(atac_bin, na.rm = TRUE)
  score_adj <- score_mul_strict / (activity_penalty + 1e-6)

  results <- rbind(results, data.frame(
    peak = peak,
    gene = gene,
    link_score = link_score,
    add = score_add,
    mul = score_mul,
    mul_weigh = score_mul_weighted,
    mul_strict = score_mul_strict,
    adj = score_adj
  ))
}

write.csv(results, "test_scores.csv", row.names = FALSE)

head(results[order(-results$link_score), ], 20)
head(results[order(-results$add), ], 20)
head(results[order(-results$mul), ], 20)
head(results[order(-results$mul_weigh), ], 20)
head(results[order(-results$mul_strict), ], 20)
head(results[order(-results$adj), ], 20)

cor(results$link_score, results$mul_weigh, use = "complete.obs")
cor(results$link_score, results$mul_strict, use = "complete.obs")
cor(results$link_score, results$adj, use = "complete.obs")
cor(results$add, results$mul_weigh, use = "complete.obs")
cor(results$add, results$mul_strict, use = "complete.obs")
cor(results$mul_strict, results$adj, use = "complete.obs")

top_link <- head(results[order(-results$link_score), c("peak", "gene")], 50)
top_add <- head(results[order(-results$add), c("peak", "gene")], 50)
top_mulw <- head(results[order(-results$mul_weigh), c("peak", "gene")], 50)
top_muls <- head(results[order(-results$mul_strict), c("peak", "gene")], 50)
top_adj <- head(results[order(-results$adj), c("peak", "gene")], 50)

nrow(merge(top_link, top_add))
nrow(merge(top_link, top_mulw))
nrow(merge(top_link, top_muls))
nrow(merge(top_link, top_adj))

## =========================
## DISTANCE-BASED EXTENSION
## =========================

## ---- 1. Build gene-coordinate table safely from annotations ----
gene_names <- mcols(annotations)$gene_name
seqs <- as.character(seqnames(annotations))
starts <- start(annotations)
ends <- end(annotations)

keep <- !is.na(gene_names) & gene_names != ""

ann_df <- data.frame(
  gene = gene_names[keep],
  gene_chr = seqs[keep],
  gene_start = starts[keep],
  gene_end = ends[keep],
  stringsAsFactors = FALSE
)

# collapse transcript/exon rows to one span per gene
gene_coord_df <- aggregate(
  cbind(gene_start, gene_end) ~ gene + gene_chr,
  data = ann_df,
  FUN = function(x) c(min = min(x), max = max(x))
)

gene_coord_df <- data.frame(
  gene = gene_coord_df$gene,
  gene_chr = gene_coord_df$gene_chr,
  gene_start = gene_coord_df$gene_start[, "min"],
  gene_end = gene_coord_df$gene_end[, "max"],
  stringsAsFactors = FALSE
)

gene_coord_df$gene_mid <- (gene_coord_df$gene_start + gene_coord_df$gene_end) / 2

## ---- 2. Parse peak coordinates from results$peak ----
peak_parts <- strsplit(results$peak, "-", fixed = TRUE)

results$peak_chr <- vapply(peak_parts, `[`, character(1), 1)
results$peak_start <- as.numeric(vapply(peak_parts, `[`, character(1), 2))
results$peak_end <- as.numeric(vapply(peak_parts, `[`, character(1), 3))
results$peak_mid <- (results$peak_start + results$peak_end) / 2

## ---- 3. Join gene coordinates ----
results <- merge(
  results,
  gene_coord_df[, c("gene", "gene_chr", "gene_start", "gene_end", "gene_mid")],
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)

## ---- 4. Compute genomic distance ----
results$distance_bp <- ifelse(
  is.na(results$gene_chr) | results$peak_chr != results$gene_chr,
  Inf,
  abs(results$peak_mid - results$gene_mid)
)

## ---- 5. Long-tail distance prior ----
d0 <- 50000

results$distance_score <- ifelse(
  is.finite(results$distance_bp),
  1 / (1 + (results$distance_bp / d0)^2),
  0
)

## ---- 6. Combine strict score + distance ----
results$final_w01 <- 0.9 * results$mul_strict + 0.1 * results$distance_score
results$final_w02 <- 0.8 * results$mul_strict + 0.2 * results$distance_score
results$final_w04 <- 0.6 * results$mul_strict + 0.4 * results$distance_score

## ---- 7. Save ----
write.csv(results, "test_scores_with_distance.csv", row.names = FALSE)

## ---- 8. Inspect ----
head(results[order(-results$final_w01), ], 20)
head(results[order(-results$final_w02), ], 20)
head(results[order(-results$final_w04), ], 20)

cor(results$link_score, results$final_w01, use = "complete.obs")
cor(results$link_score, results$final_w02, use = "complete.obs")
cor(results$link_score, results$final_w04, use = "complete.obs")

top_link <- head(results[order(-results$link_score), c("peak", "gene")], 50)
top_muls <- head(results[order(-results$mul_strict), c("peak", "gene")], 50)

top_final_w01 <- head(results[order(-results$final_w01), c("peak", "gene")], 50)
top_final_w02 <- head(results[order(-results$final_w02), c("peak", "gene")], 50)
top_final_w04 <- head(results[order(-results$final_w04), c("peak", "gene")], 50)

nrow(merge(top_link, top_final_w01))
nrow(merge(top_link, top_final_w02))
nrow(merge(top_link, top_final_w04))

nrow(merge(top_muls, top_final_w01))
nrow(merge(top_muls, top_final_w02))
nrow(merge(top_muls, top_final_w04))

head(results[order(-results$final_w02), ], 20)
cor(results$link_score, results$final_w02, use = "complete.obs")
nrow(merge(top_link, top_final_w02))
nrow(merge(top_muls, top_final_w02))

plot(
  results$link_score,
  results$final_w02,
  pch = 16,
  cex = 0.5,
  xlab = "Baseline LinkPeaks score",
  ylab = "final_w02 score",
  main = "Baseline vs final_w02"
)

results$rank_link <- rank(-results$link_score, ties.method = "average")
results$rank_final_w02 <- rank(-results$final_w02, ties.method = "average")
results$rank_diff <- results$rank_link - results$rank_final_w02

head(results[order(results$rank_diff), c("peak", "gene", "link_score", "final_w02", "rank_link", "rank_final_w02", "rank_diff")], 20)
head(results[order(-results$rank_diff), c("peak", "gene", "link_score", "final_w02", "rank_link", "rank_final_w02", "rank_diff")], 20)

top_link_50 <- head(results[order(-results$link_score), ], 50)
top_final_50 <- head(results[order(-results$final_w02), ], 50)

summary(top_link_50$distance_bp)
summary(top_final_50$distance_bp)

median(top_link_50$distance_bp, na.rm = TRUE)
median(top_final_50$distance_bp, na.rm = TRUE)

top_link_20 <- head(results[order(-results$link_score), c("peak", "gene", "link_score", "mul_strict", "distance_bp", "distance_score", "final_w02")], 20)
top_final_20 <- head(results[order(-results$final_w02), c("peak", "gene", "link_score", "mul_strict", "distance_bp", "distance_score", "final_w02")], 20)

top_link_20
top_final_20

hist(results$distance_bp, breaks = 50)
hist(results$distance_bp[order(-results$final_w02)[1:200]], breaks = 50)

mean(top_final_50$distance_bp > 50000)

results$final_test <- results$mul_strict * (0.5 + 0.5 * results$distance_score)

top_test_50 <- head(results[order(-results$final_test), ], 50)
mean(top_test_50$distance_bp > 50000)
median(top_test_50$distance_bp)

results$final_v4 <- results$mul_weigh * (0.5 + 0.5 * results$distance_score)

top_v4_50 <- head(results[order(-results$final_v4), ], 50)

mean(top_v4_50$distance_bp > 50000)
median(top_v4_50$distance_bp)

results$final_v5 <- results$mul_weigh * (0.7 + 0.3 * results$distance_score)

top_v5_50 <- head(results[order(-results$final_v5), ], 50)
mean(top_v5_50$distance_bp > 50000)
median(top_v5_50$distance_bp)

top_v5_50[top_v5_50$distance_bp > 50000, 
          c("gene", "peak", "distance_bp", "mul_weigh", "final_v5")]


top_link_50 <- head(results[order(-results$link_score), c("peak","gene")], 50)
top_v5_50 <- head(results[order(-results$final_v5), c("peak","gene")], 50)

nrow(merge(top_link_50, top_v5_50))

for (n in c(10, 20, 50, 100, 200)) {
  top_link_n <- head(results[order(-results$link_score), c("peak","gene")], n)
  top_v5_n <- head(results[order(-results$final_v5), c("peak","gene")], n)
  cat("Top", n, "overlap:", nrow(merge(top_link_n, top_v5_n)), "\n")
}

cor(results$link_score, results$final_v5, use = "complete.obs")

results$rank_v5 <- rank(-results$final_v5)

results$rank_diff_v5 <- results$rank_link - results$rank_v5

# promoted by your method
head(results[order(results$rank_diff_v5), 
     c("gene","peak","link_score","final_v5","distance_bp","rank_link","rank_v5")], 20)

# demoted by your method
head(results[order(-results$rank_diff_v5), 
     c("gene","peak","link_score","final_v5","distance_bp","rank_link","rank_v5")], 20)

write.csv(results, "test_scores_with_distance_FINAL.csv", row.names = FALSE)


######################

## =========================
## LOAD SAVED RESULTS (NO MEMORY DEPENDENCE)
## =========================

# set path if needed
file_path <- "test_scores_with_distance_FINAL.csv"

# check file exists
if (!file.exists(file_path)) {
  stop("ERROR: results file not found. Check path.")
}

# read results
results <- read.csv(file_path, stringsAsFactors = FALSE)

cat("Loaded results with", nrow(results), "rows\n")

# verify required columns exist
required_cols <- c("gene", "link_score", "final_v5")

missing_cols <- setdiff(required_cols, colnames(results))
if (length(missing_cols) > 0) {
  stop(paste("ERROR: Missing columns:", paste(missing_cols, collapse = ", ")))
}

cat("Columns OK. final_v5 present.\n")

## optional: quick sanity check
summary(results$final_v5)
summary(results$link_score)

######ORA TEST

## =========================
## DEFINE GENE SETS
## =========================


## ---- 1. choose top genes ----
top_n <- 100

baseline_genes <- unique(head(results[order(-results$link_score), "gene"], top_n))
v5_genes <- unique(head(results[order(-results$final_v5), "gene"], top_n))

# background = all genes present in your candidate table
background_genes <- unique(results$gene)

cat("Baseline top genes:", length(baseline_genes), "\n")
cat("V5 top genes:", length(v5_genes), "\n")
cat("Background genes:", length(background_genes), "\n")

## ---- 2. map gene symbols to Entrez IDs ----
baseline_map <- bitr(
  baseline_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

v5_map <- bitr(
  v5_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

background_map <- bitr(
  background_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

baseline_entrez <- unique(baseline_map$ENTREZID)
v5_entrez <- unique(v5_map$ENTREZID)
background_entrez <- unique(background_map$ENTREZID)

cat("Baseline mapped Entrez IDs:", length(baseline_entrez), "\n")
cat("V5 mapped Entrez IDs:", length(v5_entrez), "\n")
cat("Background mapped Entrez IDs:", length(background_entrez), "\n")

## ---- 3. run ORA (GO BP) ----
ora_baseline <- enrichGO(
  gene = baseline_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

ora_v5 <- enrichGO(
  gene = v5_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

## ---- 4. save ORA tables ----
baseline_ora_df <- as.data.frame(ora_baseline)
v5_ora_df <- as.data.frame(ora_v5)

write.csv(baseline_ora_df, "ora_baseline_GO_BP.csv", row.names = FALSE)
write.csv(v5_ora_df, "ora_final_v5_GO_BP.csv", row.names = FALSE)

## ---- 5. inspect top terms in console ----
cat("\nTop baseline GO terms:\n")
print(head(baseline_ora_df[, c("Description", "GeneRatio", "BgRatio", "p.adjust", "Count")], 15))

cat("\nTop final_v5 GO terms:\n")
print(head(v5_ora_df[, c("Description", "GeneRatio", "BgRatio", "p.adjust", "Count")], 15))

## ---- 6. make separate dotplots ----
p_baseline <- dotplot(ora_baseline, showCategory = 15) +
  ggtitle("ORA: Baseline (top genes)") +
  theme_bw()

p_v5 <- dotplot(ora_v5, showCategory = 15) +
  ggtitle("ORA: final_v5 (top genes)") +
  theme_bw()

print(p_baseline)
print(p_v5)

ggsave("ora_baseline_dotplot.png", p_baseline, width = 10, height = 7, dpi = 300)
ggsave("ora_final_v5_dotplot.png", p_v5, width = 10, height = 7, dpi = 300)

## ---- 7. simple side-by-side comparison plot ----
baseline_top <- head(baseline_ora_df, 10)
baseline_top$method <- "baseline"

v5_top <- head(v5_ora_df, 10)
v5_top$method <- "final_v5"

ora_compare <- rbind(
  baseline_top[, c("Description", "Count", "p.adjust", "method")],
  v5_top[, c("Description", "Count", "p.adjust", "method")]
)

ora_compare$log10_padj <- -log10(ora_compare$p.adjust)

# keep term order readable
ora_compare$Description <- factor(
  ora_compare$Description,
  levels = rev(unique(ora_compare$Description))
)

p_compare <- ggplot(
  ora_compare,
  aes(x = Count, y = Description, color = method, size = log10_padj)
) +
  geom_point() +
  facet_wrap(~method, scales = "free_y") +
  theme_bw() +
  ggtitle("GO BP ORA: baseline vs final_v5")

print(p_compare)

ggsave("ora_baseline_vs_final_v5.png", p_compare, width = 12, height = 8, dpi = 300)

nrow(as.data.frame(ora_baseline))
nrow(as.data.frame(ora_v5))

length(baseline_genes)
length(v5_genes)

length(unique(baseline_genes))
length(unique(v5_genes))


head(baseline_map)
head(v5_map)

length(baseline_entrez)
length(v5_entrez)


for (n in c(50, 100, 200)) {
  cat("\nTop", n, "\n")
  
  baseline_genes <- unique(head(results[order(-results$link_score), "gene"], n))
  v5_genes <- unique(head(results[order(-results$final_v5), "gene"], n))
  
  cat("Baseline unique:", length(baseline_genes), "\n")
  cat("V5 unique:", length(v5_genes), "\n")
}
