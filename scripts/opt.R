opt <- list(
  data_dir = "data",
  h5_file = "filtered_feature_bc_matrix.h5",
  frag_file = "atac_fragments.tsv.gz",

  candidate_top_k = 10000L,
  link_distance = 500000L,
  distance_d0 = 50000,
  lambda_distance = 0.3, #0.3
  alpha_tf = 0.5, # 0.5

  cluster_resolution = 0.5,
  pca_dims = 30L,
  lsi_dims_start = 2L,
  lsi_dims_end = 30L,

  species = 9606L,
  collection = "CORE",
  motif_min_score = NULL,
  tf_expressed_frac = 0.10,
  top_k_motif_names = 3L,

  ora_top_n = 100L,
  ora_show_category = 15L,

  tier_high_quantile = 0.90,
  tier_medium_quantile = 0.70,

  output_prefix = file.path("results", "multiome_rie")
)