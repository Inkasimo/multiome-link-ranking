# Split reranker workflow

This is a proposed split of the original monolithic reranker script into two scripts:

```text
scripts/run_linkpeaks_reranker.R
scripts/evaluate_rankings.R
```

The heavy script generates candidate links and feature columns once. The light script consumes that feature table to produce default rankings, ablations, ORA tables, plots, and metrics.

## Run locally

```bash
snakemake --snakefile workflow/Snakefile --configfile config/default.yaml --cores 4
```

## Run inside Docker

```bash
docker run --rm \
  -v "$PWD":/work \
  -w /work \
  multiome-link-ranking:pilot \
  snakemake --snakefile workflow/Snakefile --configfile config/default.yaml --cores 4
```

## Heavy step only

```bash
Rscript scripts/run_linkpeaks_reranker.R \
  --data-dir data \
  --h5-file filtered_feature_bc_matrix.h5 \
  --frag-file atac_fragments.tsv.gz \
  --candidate-top-k 5000 \
  --link-distance 500000 \
  --distance-d0 50000 \
  --output-dir results/pbmc/features \
  --run-name pbmc
```

## Evaluate full model only

```bash
Rscript scripts/evaluate_rankings.R \
  --features-file results/pbmc/features/pbmc_link_features.csv \
  --baseline-file results/pbmc/features/pbmc_baseline_links_full.csv \
  --baseline-distance-file results/pbmc/features/pbmc_baseline_links_with_distance.csv \
  --score-mode full \
  --lambda-distance 0.30 \
  --alpha-tf 0.50 \
  --output-dir results/pbmc/rankings/full \
  --run-name pbmc_full
```

## Scoring modes

```text
linkpeaks
coactivity
distance_only
coactivity_distance
coactivity_tf
full
```
