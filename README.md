# Multiome Link Ranking

Interpretable reranking of candidate peak–gene links from single-cell multiome data.

## Overview

This repository contains a prototype workflow for ranking candidate regulatory peak–gene links from paired single-cell RNA + ATAC multiome data.

The current workflow is a **reranker benchmark**. It uses `Signac::LinkPeaks()` to generate an initial candidate set of peak–gene links, then reranks those candidates using a biologically structured score based on:

1. RNA–ATAC coactivity
2. genomic distance to the gene transcription start site
3. transcription factor / motif support

The goal of this branch is not yet to provide a fully standalone peak–gene linking model. Instead, it tests whether the scoring idea contains useful biological signal before replacing LinkPeaks candidate generation in a later standalone version.

In short:

```text
current branch:
  LinkPeaks candidates -> interpretable reranking -> benchmark outputs

future standalone model:
  own cis-window candidates -> same scoring logic -> independent peak–gene linker
```

## Biological idea

A regulatory peak–gene link should be prioritized when three signals agree:

* the peak is accessible in cells where the gene is expressed
* the peak is genomically plausible for the gene
* the peak contains motif evidence for transcription factors that are expressed

The score is deliberately interpretable. It is not a black-box model and does not learn unconstrained weights. Instead, it combines a small number of biologically motivated terms.

The core assumption is:

```text
regulatory support = coactivity × distance prior × TF/motif support
```

This makes the method useful as a transparent benchmarkable prototype before moving to a full standalone cis-window model.

## Current workflow

The main R workflow:

1. reads 10x multiome RNA + ATAC input
2. builds a Seurat object
3. creates a Signac chromatin assay
4. performs RNA preprocessing
5. performs ATAC preprocessing
6. integrates modalities using WNN
7. clusters cells
8. runs `LinkPeaks()` to generate baseline candidate peak–gene links
9. computes reranking scores
10. adds distance and motif/TF support
11. writes ranked links and benchmark summaries

The current reranker starts from LinkPeaks candidates. Therefore, it should be interpreted as a **candidate reranking method**, not yet a de novo candidate-generation method.

## Scoring formula

For a candidate peak–gene pair ((p,g)), let:

* (C) be the set of cells
* (x_{gc}^{RNA}) be normalized RNA expression of gene (g) in cell (c)
* (x_{pc}^{ATAC}) be normalized ATAC accessibility of peak (p) in cell (c)
* (d_{pg}) be the distance between the peak midpoint and gene TSS
* (T_p) be the TF/motif support score for peak (p)
* (d_0) be the distance decay scale
* (\lambda) be the distance strength
* (\alpha) be the TF/motif strength

### 1. Standardized RNA and ATAC signals

For each gene and peak, values are z-scored across cells:

```text
z_gc^RNA  = (x_gc^RNA  - mean(x_g^RNA))  / sd(x_g^RNA)
z_pc^ATAC = (x_pc^ATAC - mean(x_p^ATAC)) / sd(x_p^ATAC)
```

If the standard deviation is zero, the z-score vector is set to zero.

### 2. Positive RNA–ATAC coactivity

The core coactivity term is:

```text
A_pg = (1 / |C|) * Σ_c max(z_gc^RNA, 0) * max(z_pc^ATAC, 0)
```

This is the `mul_weigh` score.

Only cells with positive RNA and positive ATAC deviations contribute. This makes the term behave like a soft biological AND gate: a peak–gene pair is rewarded when accessibility and expression are high together in the same cells.

### 3. Distance prior

The distance prior is a smooth long-tail decay:

```text
D_pg = 1 / (1 + (d_pg / d0)^2)
```

Nearby links receive values close to 1. Distant links are downweighted but not removed completely.

The distance-regularized coactivity score is:

```text
final_v5 = A_pg * ((1 - λ) + λ * D_pg)
```

This means distance acts as a prior, not as a hard cutoff.

### 4. Continuous TF/motif support

For each peak, motif support is computed using continuous motif scores and expression of the corresponding transcription factors.

Conceptually:

```text
T_p = scaled Σ_m motif_score_pm * TF_expression_m
```

where:

* `motif_score_pm` is the continuous motif score for motif (m) in peak (p)
* `TF_expression_m` is the expression-based weight for the transcription factor associated with motif (m)

This gives higher support to peaks containing strong motifs for TFs that are actually expressed.

### 5. Final reranking score

The final score is:

```text
S_pg = A_pg * ((1 - λ) + λ * D_pg) * (1 + α * T_p)
```

Equivalently:

```text
final_v6 = mul_weigh *
           ((1 - lambda_distance) + lambda_distance * distance_score) *
           (1 + alpha_tf * tf_score)
```

Higher (S_{pg}) means a higher-priority candidate regulatory link.

## Main interpretation

The three score components have distinct roles:

```text
coactivity:
  main RNA–ATAC evidence

distance:
  smooth genomic plausibility prior

TF/motif:
  regulatory plausibility support
```

The score is multiplicative because weak evidence in one component should reduce confidence, while strong agreement across components should promote the link.

## Current benchmark position

This workflow currently compares or supports comparison against:

* LinkPeaks baseline
* coactivity-only ranking
* distance-only baseline
* full reranker
* SCENT where feasible
* ArchR attempted but excluded from the prototype benchmark if no usable links are produced

The important baseline checks are:

```text
Does the full reranker beat LinkPeaks?
Does it beat coactivity alone?
Does it avoid collapsing into distance-only behavior?
Does TF/motif support improve biological coherence?
Does it preserve plausible distal links?
```

## Expected outputs

Typical output files include:

```text
results/*_baseline_links_full.csv
results/*_test_scores.csv
results/*_ranked_links.csv
results/*_top100_final_links.csv
results/*_top_promoted_links.csv
results/*_top_demoted_links.csv
results/*_tier_summary.csv
results/*_summary_metrics.csv
results/*_ora_final_v6_GO_BP.csv
results/*_distance_only_summary.csv
results/*_distance_only_top100_links.csv
results/*_distance_only_ora.csv
results/*_distance_distribution.png
results/*_baseline_vs_final_v6_scatter.png
```

The main output is:

```text
results/*_ranked_links.csv
```

with columns such as:

```text
peak
gene
link_score
mul_weigh
distance_bp
distance_score
tf_score
final_v6
tier
rank_final_v6
```

## Docker usage

Build the Docker image:

```bash
docker build -t multiome-link-ranking:pilot -f containers/Dockerfile .
```

Smoke-test the R stack:

```bash
docker run --rm multiome-link-ranking:pilot Rscript -e 'library(Seurat); library(Signac); library(TFBSTools); library(JASPAR2022); library(motifmatchr); cat("R stack OK\n")'
```

Run from the repository root by mounting the project into the container:

```bash
docker run --rm \
  -v "$PWD":/work \
  -w /work \
  multiome-link-ranking:pilot \
  Rscript scripts/multiome_rie_v1.R \
    --data-dir data \
    --h5-file filtered_feature_bc_matrix.h5 \
    --frag-file atac_fragments.tsv.gz \
    --candidate-top-k 5000 \
    --link-distance 500000 \
    --distance-d0 50000 \
    --lambda-distance 0.30 \
    --alpha-tf 0.50 \
    --output-prefix results/multiome_rie
```

The fragment index is required:

```text
atac_fragments.tsv.gz.tbi
```

## Current limitations

This branch is a proof-of-concept reranker, not the final standalone method.

Important limitations:

* candidate links are still generated by LinkPeaks
* the method does not yet generate its own full cis-window candidate universe
* current scoring is global rather than cluster/metacell-aware
* external validation is still limited
* ArchR integration is deferred unless a stable benchmarkable output is produced
* SCENT is computationally heavy and should be run on tractable subsets or chromosome-level benchmarks

## Planned standalone extension

The next methodological step is to replace LinkPeaks candidate generation with an internal cis-window candidate builder.

Planned standalone candidate generation:

```text
for each gene:
  keep peaks on the same chromosome
  keep peaks within a fixed cis window
  apply gene expression prevalence filter
  apply peak accessibility prevalence filter
```

Initial candidate windows:

```text
100–500 kb
```

Later sensitivity:

```text
1 Mb
```

The first standalone version should keep the current scoring formula unchanged:

```text
score = coactivity × distance prior × TF/motif support
```

This tests the key question:

```text
Does the scoring model still work when LinkPeaks no longer defines the candidate universe?
```

## Future cluster/metacell-aware model

After the standalone cis-window version works, the model can be extended to cluster-aware or metacell-based scoring.

For cluster (k):

```text
S_pgk = A_pgk × D_pg × (1 + α T_pgk)
```

where:

* (A_{pgk}) is coactivity within cluster/metacells
* (D_{pg}) is the global distance prior
* (T_{pgk}) is cluster-aware TF/motif support

Then aggregate across clusters:

```text
S_pg = max_k S_pgk
best_cluster = argmax_k S_pgk
```

This would allow the model to identify regulatory links that are active only in specific cell states.

## Project status

Current branch:

```text
workflow-benchmark-pilot
```

Purpose:

```text
Dockerized, reproducible reranker benchmark
```

Ultimate project goal:

```text
standalone interpretable peak–gene linking model
```

The current benchmark branch is the evidence-gathering step. It tests whether the scoring idea is strong enough to justify building the standalone model.
