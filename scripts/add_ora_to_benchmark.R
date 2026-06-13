#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
})

option_list <- list(
  make_option(c("--benchmark-dir"), type = "character",
              default = file.path("results", "benchmark_reranker_scent_sweep"),
              help = "Directory containing benchmark_all_methods_combined.csv or benchmark_*_links.csv [default %default]"),
  make_option(c("--top-n"), type = "integer", default = 100,
              help = "Top N links/genes per method for GO ORA [default %default]"),
  make_option(c("--ora-show-category"), type = "integer", default = 15,
              help = "Number of GO terms shown in dotplots [default %default]"),
  make_option(c("--pvalue-cutoff"), type = "double", default = 0.05,
              help = "enrichGO pvalueCutoff [default %default]"),
  make_option(c("--qvalue-cutoff"), type = "double", default = 0.20,
              help = "enrichGO qvalueCutoff [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
names(opt) <- gsub("-", "_", names(opt), fixed = TRUE)

msg <- function(...) cat(sprintf(...), "\n")

if (!dir.exists(opt$benchmark_dir)) {
  stop("Benchmark directory not found: ", opt$benchmark_dir)
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
  entrez_ids <- unique(entrez_ids)
  universe_ids <- unique(universe_ids)

  if (length(entrez_ids) == 0 || length(universe_ids) == 0) {
    return(NULL)
  }

  tryCatch({
    enrichGO(
      gene = entrez_ids,
      universe = universe_ids,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = opt$pvalue_cutoff,
      qvalueCutoff = opt$qvalue_cutoff,
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
    xlim(0, 2) +
    ylim(0, 2) +
    theme_void() +
    ggtitle(title_txt)
}

make_dotplot_safe <- function(ora_obj, title_txt, show_n = 15) {
  if (is.null(ora_obj) || nrow(as.data.frame(ora_obj)) == 0) {
    return(empty_plot(title_txt))
  }

  dotplot(ora_obj, showCategory = show_n) +
    ggtitle(title_txt) +
    theme_bw()
}

read_benchmark_links <- function(benchmark_dir) {
  combined_path <- file.path(benchmark_dir, "benchmark_all_methods_combined.csv")

  if (file.exists(combined_path)) {
    msg("Reading combined benchmark table: %s", combined_path)
    x <- fread(combined_path)

    required <- c("peak", "gene", "score", "method")
    missing <- setdiff(required, names(x))
    if (length(missing) > 0) {
      stop("Combined benchmark table missing columns: ", paste(missing, collapse = ", "))
    }

    return(x[, .(peak, gene, score, method)])
  }

  msg("Combined table not found. Reading individual benchmark_*_links.csv files.")

  files <- list.files(
    benchmark_dir,
    pattern = "^benchmark_.*_links\\.csv$",
    full.names = TRUE
  )

  files <- files[!grepl("all_methods_combined", basename(files))]
  if (length(files) == 0) {
    stop("No benchmark link tables found in: ", benchmark_dir)
  }

  rbindlist(lapply(files, function(f) {
    x <- fread(f)
    required <- c("peak", "gene", "score", "method")
    missing <- setdiff(required, names(x))
    if (length(missing) > 0) {
      stop("File missing columns: ", f, " / ", paste(missing, collapse = ", "))
    }
    x[, .(peak, gene, score, method)]
  }), use.names = TRUE, fill = TRUE)
}

links <- read_benchmark_links(opt$benchmark_dir)
links <- links[!is.na(gene) & nzchar(gene) & is.finite(score)]

if (nrow(links) == 0) {
  stop("No usable benchmark links found after filtering.")
}

methods <- sort(unique(links$method))
msg("Methods found: %s", paste(methods, collapse = ", "))

background_genes <- unique(links$gene)
background_map <- safe_bitr(background_genes)
background_entrez <- unique(background_map$ENTREZID)

msg("Background genes: %d", length(background_genes))
msg("Background mapped Entrez IDs: %d", length(background_entrez))

ora_rows <- list()
top_gene_rows <- list()

for (method_name in methods) {
  msg("Running GO BP ORA for %s", method_name)

  method_dt <- links[method == method_name]
  setorder(method_dt, -score)

  top_genes <- unique(method_dt[1:min(.N, opt$top_n)]$gene)
  top_map <- safe_bitr(top_genes)
  top_entrez <- unique(top_map$ENTREZID)

  top_gene_rows[[method_name]] <- data.table(
    method = method_name,
    rank = seq_along(top_genes),
    gene = top_genes
  )

  ora <- safe_enrich_go(top_entrez, background_entrez)
  ora_df <- if (is.null(ora)) data.frame() else as.data.frame(ora)

  fwrite(
    as.data.table(ora_df),
    file.path(opt$benchmark_dir, sprintf("benchmark_%s_ora_GO_BP.csv", method_name))
  )

  p <- make_dotplot_safe(
    ora,
    sprintf("GO BP ORA: %s top %d genes", method_name, opt$top_n),
    opt$ora_show_category
  )

  ggsave(
    file.path(opt$benchmark_dir, sprintf("benchmark_%s_ora_dotplot.png", method_name)),
    p,
    width = 10,
    height = 7,
    dpi = 300
  )

  ora_rows[[method_name]] <- data.table(
    method = method_name,
    top_n = opt$top_n,
    top_genes = length(top_genes),
    mapped_top_entrez = length(top_entrez),
    background_genes = length(background_genes),
    mapped_background_entrez = length(background_entrez),
    n_ora_terms = nrow(ora_df)
  )
}

ora_summary <- rbindlist(ora_rows, use.names = TRUE, fill = TRUE)
top_gene_table <- rbindlist(top_gene_rows, use.names = TRUE, fill = TRUE)

fwrite(ora_summary, file.path(opt$benchmark_dir, "benchmark_ora_summary.csv"))
fwrite(top_gene_table, file.path(opt$benchmark_dir, sprintf("benchmark_top%d_genes_for_ORA.csv", opt$top_n)))

p_count <- ggplot(ora_summary, aes(x = reorder(method, n_ora_terms), y = n_ora_terms)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    title = sprintf("GO BP ORA terms by method: top %d genes", opt$top_n),
    x = NULL,
    y = "Significant GO BP terms"
  )

ggsave(
  file.path(opt$benchmark_dir, sprintf("benchmark_ora_term_counts_top%d.png", opt$top_n)),
  p_count,
  width = 9,
  height = 6,
  dpi = 300
)

# Merge ORA counts into benchmark_summary_metrics.csv if present.
summary_path <- file.path(opt$benchmark_dir, "benchmark_summary_metrics.csv")
if (file.exists(summary_path)) {
  summary_dt <- fread(summary_path)
  summary_dt[, n_ora_terms := NULL]
  summary_dt <- merge(summary_dt, ora_summary[, .(method, n_ora_terms)], by = "method", all.x = TRUE)
  fwrite(summary_dt, summary_path)
  msg("Updated benchmark_summary_metrics.csv with n_ora_terms.")
}

msg("ORA complete.")
msg("ORA summary: %s", file.path(opt$benchmark_dir, "benchmark_ora_summary.csv"))
