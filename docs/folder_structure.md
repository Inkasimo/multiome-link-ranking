# Proposed folder structure

```text
config/
  default.yaml              # dataset paths and default parameters
  ablations.yaml            # scoring modes and ablation parameters

scripts/
  run_linkpeaks_reranker.R  # heavy: object build, LinkPeaks, coactivity, distance, motif/TF features
  evaluate_rankings.R       # light: default ranking, ablations, ORA, plots, metrics

workflow/
  Snakefile                 # one heavy feature rule, many lightweight evaluation rules

results/
  <dataset>/
    features/
      <dataset>_link_features.csv
      <dataset>_baseline_links_full.csv
      <dataset>_baseline_links_with_distance.csv
      .done

    rankings/
      linkpeaks/
      coactivity/
      coactivity_distance/
      coactivity_tf/
      distance_only/
      full/
```

## Design rule

Run the heavy Seurat/Signac/LinkPeaks/motif workflow once. It writes a canonical feature table.

Run every ablation from that feature table. This avoids rebuilding the object, rerunning LinkPeaks, and recomputing motif scores for every score variant.
