opt <- list(
  results_dir = file.path("results", "alpha_tf_05_after_tss_current"),
  baseline_file = "multiome_rie_baseline_links_full.csv",
  ranked_file = "multiome_rie_ranked_links.csv",
  test_scores_file = "multiome_rie_test_scores.csv",

  data_dir = "data",
  h5_file = "filtered_feature_bc_matrix.h5",
  frag_file = "atac_fragments.tsv.gz",

  output_dir = file.path("results", "benchmark_panel"),
  sample_name = "multiome_sample",
  genome = "hg38",

  link_distance = 500000,
  distance_d0 = 50000,
  distal_threshold = 50000,

  top_n = 100,
  top_k_compare = 200,

  min_pair_frac = 0.05,
  scoring_celltype = NA_character_,

  archr_threads = 4L,
  scent_cores = 4L,
  scent_regr = "poisson",

  ora_show_category = 15L,
  seed = 42L,

  skip_archr = FALSE,
  skip_scent = FALSE,
  archr_force = FALSE
)

set.seed(opt$seed)
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)