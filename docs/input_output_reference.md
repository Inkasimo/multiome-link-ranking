# Input / output reference

Practical interface documentation for the PBMC LinkPeaks-candidate peak–gene reranking
benchmark. Column semantics are read from the committed output files; parameter semantics
are read from the four pipeline scripts and the Snakefile.

Paths use `<dataset>` = the `dataset` key in `config/default.yaml`, which is `pbmc` for all
committed results.

---

## 1. Required input files

All three live under `data/` (`data_dir` in `config/default.yaml`) and are **not** committed
to the repository.

| File | Config key | Purpose | Consumed by |
|---|---|---|---|
| `data/filtered_feature_bc_matrix.h5` | `h5_file` | 10x multiome filtered feature-barcode matrix, containing both the Gene Expression and Peaks assays | `run_linkpeaks_reranker.R`, `run_scent_chr_sweep.R` |
| `data/atac_fragments.tsv.gz` | `frag_file` | ATAC fragment file | `run_linkpeaks_reranker.R` |
| `data/atac_fragments.tsv.gz.tbi` | derived: `frag_file + ".tbi"` | Tabix index. Declared as an explicit Snakemake input, so a missing index fails at DAG construction rather than mid-run. | `run_linkpeaks_reranker.R` |

### Additional required resource

| File | Purpose |
|---|---|
| `resources/jaspar/JASPAR2022.sqlite` | Local JASPAR2022 motif database |
| `resources/jaspar/JASPAR2022.sqlite.sha256` | Checksum |

This is a hard requirement, not an optimisation. `scripts/run_linkpeaks_reranker.R` lines
8–48 define `seed_jaspar2022_cache()`, which registers the local file with `BiocFileCache`
under the remote URL `https://jaspar2022.genereg.net/download/database/JASPAR2022.sqlite`.
Without it the `JASPAR2022` package attempts a network download at motif-loading time. Because
the container is intended to run offline against a bind-mounted working directory, this
download fails and feature generation aborts. Pinning this file is what unblocked the
recorded feature-generation run.

The file must be present in the **bind-mounted working directory**, not only in the image —
the seeding function resolves it via `file.path(getwd(), "resources", "jaspar", ...)`.

### Reference annotation

Not files, but package dependencies resolved from `renv.lock`:

- `EnsDb.Hsapiens.v86` — gene and transcript annotation, TSS derivation
- `BSgenome.Hsapiens.UCSC.hg38` — sequence for `RegionStats()` and motif matching

Genome build is **hg38**, set explicitly at `run_linkpeaks_reranker.R` line 377
(`genome(annotations) <- "hg38"`) and line 384.

---

## 2. Required config fields

### `config/default.yaml`

| Field | Committed value | Meaning |
|---|---|---|
| `dataset` | `pbmc` | Run name; determines output paths and output-file prefixes |
| `data_dir` | `data` | Input directory |
| `h5_file` | `filtered_feature_bc_matrix.h5` | Filename within `data_dir` |
| `frag_file` | `atac_fragments.tsv.gz` | Filename within `data_dir`; `.tbi` inferred |
| `output_root` | `results` | Root of the output tree |
| `default_score_mode` | `full` | Mode used by `run_default_score` |
| `candidate_top_k` | `5000` | Number of LinkPeaks candidates retained for the feature table |
| `link_distance` | `500000` | Maximum peak-to-TSS distance passed to `LinkPeaks()` |
| `distance_d0` | `50000` | Decay scale $d_0$ in the distance score |
| `cluster_resolution` | `0.5` | Louvain resolution on the WNN graph |
| `pca_dims` | `30` | RNA PCA dimensions |
| `lsi_dims_start` | `2` | First LSI component used (component 1 dropped — depth-correlated) |
| `lsi_dims_end` | `30` | Last LSI component used |
| `species` | `9606` | NCBI taxonomy ID for JASPAR motif retrieval |
| `collection` | `CORE` | JASPAR collection |
| `tf_expressed_frac` | `0.10` | Minimum fraction of cells expressing a TF for its motif to be retained |
| `top_k_motif_names` | `3` | Number of motif names recorded per peak in `peak_top_motifs` |
| `seed` | `42` | RNG seed, passed to all four scripts |
| `save_checkpoints` | `false` | Controls the `--save-checkpoints` flag |
| `lambda_distance` | `0.30` | Default $\lambda$; **overridden per mode** by `config/ablations.yaml` |
| `alpha_tf` | `0.50` | Default $\alpha$; overridden per mode |
| `ora_top_n` | `100` | Top-N links whose genes form the ORA foreground |
| `ora_show_category` | `15` | Categories drawn in the ORA dotplot |
| `tier_high_quantile` | `0.90` | Score quantile above which a link is tier `High` |
| `tier_medium_quantile` | `0.70` | Score quantile above which a link is tier `Medium` |
| `ablations_file` | `config/ablations.yaml` | Score-mode definitions |

**Not currently exposed:** `--candidate-filter` (default `positive_score` in
`run_linkpeaks_reranker.R` line 79). This is the filter that reduces the raw `LinkPeaks()`
output to the 15,806 retained candidates, and it is invisible in the config. It should be
promoted to `config/default.yaml`.

### `config/ablations.yaml`

Maps a mode name to a `score_mode` plus optional `lambda_distance` / `alpha_tf` overrides.
Any mode key becomes a valid `--mode` argument and a subdirectory under
`results/<dataset>/rankings/`. Eleven modes are defined; the ten committed ones are listed
in `docs/results_report.md` §2.

### `config/scent_run.yaml`

| Field | Committed value | Meaning |
|---|---|---|
| `output_dir` | `results/pbmc/scent_chr_sweep_100kb_frac020_1000cells` | Sweep output root |
| `chromosomes` | `chr1`…`chr22` | Explicit list, or `auto` |
| `link_distance` | `100000` | Cis window for candidate generation — **note this is 100 kb, not the 500 kb used for LinkPeaks candidates** |
| `min_pair_frac` | `0.02` | Minimum fraction of cells with nonzero count, applied to genes and peaks independently |
| `max_cells` | `1000` | Cell subsample; `<= 0` uses all cells |
| `max_scent_candidates` | `100000` | Per-chromosome guardrail; the run aborts above this |
| `scent_cores` | `4` | Passed to `SCENT_algorithm()` |
| `scent_regr` | `poisson` | Regression family |
| `scoring_celltype` | `""` | Empty means a synthetic `all_cells` label is used — **no real cell-type stratification was performed** |
| `skip_existing` | `true` | Reuse existing per-chromosome outputs |

### `config/scent_validation.yaml`

| Field | Committed value | Meaning |
|---|---|---|
| `scent_sweep_dir` | `results/pbmc/scent_chr_sweep_100kb_frac020_1000cells` | Source of SCENT results |
| `output_dir` | `results/pbmc/scent_validation` | Validation output root |
| `methods` | `linkpeaks, coactivity, coactivity_tf, full_moddist_lambda_0_1, full_lambda_0_1, distance_only` | Six of the eleven modes |
| `top_n_values` | `[50, 100, 200]` | Top-N cut points |
| `top_k_compare` | `200` | Depth for pairwise overlap and the `top200_*` exports |
| `reciprocal_overlap` | `0.5` | Minimum reciprocal overlap for a ranked peak to match a SCENT peak |
| `distal_threshold` | `50000` | Threshold defining `distal_frac` |
| `distance_d0` | `50000` | Recomputation scale, kept consistent with feature generation |
| `restrict_to_scent_chrs` | `true` | Restrict ranked methods to chromosomes SCENT covered. **This is why the evaluated universe is 4,976 pairs, not 5,000.** |
| `write_standardized_links` | `false` | Suppresses an additional per-method export |

Three keys are consumed by the Snakefile but **absent** from this file, surviving only via
defaults at `workflow/Snakefile` lines 122–123:

| Field | Snakefile default | Meaning |
|---|---|---|
| `scent_support_rule` | `pvalue_positive` | Support rule: `beta > 0` and `boot_basic_p <= threshold` |
| `scent_p_threshold` | `0.05` | Threshold for the above |
| `scent_min_score` | `0.0` | Used only by the `positive_score` rule |

`config/scent_validation_with_producer.yaml` defines all three, with comments.

---

## 3. Expected directory structure

```
data/                                     # not committed
  filtered_feature_bc_matrix.h5
  atac_fragments.tsv.gz
  atac_fragments.tsv.gz.tbi

resources/jaspar/                         # not committed (size / licensing)
  JASPAR2022.sqlite
  JASPAR2022.sqlite.sha256

config/
  default.yaml
  ablations.yaml
  scent_run.yaml
  scent_validation.yaml

scripts/
  run_linkpeaks_reranker.R                # heavy: object build → LinkPeaks → features
  evaluate_rankings.R                     # light: one score mode → ranking + diagnostics
  run_scent_chr_sweep.R                   # SCENT producer, per chromosome
  benchmark_scent_validation.R            # SCENT consumer, cross-method comparison
  summarize_scent_validation_min_distance.R   # proximal-removal controls

workflow/Snakefile
containers/Dockerfile
run_analysis.py

results/<dataset>/
  features/
    .done
    <dataset>_link_features.csv
    <dataset>_baseline_links_full.csv
    <dataset>_baseline_links_with_distance.csv
  rankings/
    <mode>/
      .done
      <dataset>_<mode>_ranked_links.csv
      ... per-mode diagnostics
  scent_chr_sweep_<tag>/
    .done
    scent_chr_sweep_summary.csv
    scent_links_all_chromosomes.csv
    chr<N>/
      scent_candidates_chr<N>.csv
      scent_links_chr<N>.csv
  scent_validation/
    .done
    ... cross-method summaries
  scent_validation_min_distance/
    ... proximal-removal summaries
```

`.done` files are empty Snakemake sentinels. The workflow depends on them; do not delete
them from a working tree.

---

## 4. Important intermediate outputs

### `results/<dataset>/features/`

Written once by `run_linkpeaks_reranker.R`. Everything downstream reads from here, which is
the central design decision: the heavy Seurat/Signac/motif work runs once, and all ten score
modes are evaluated from the same table.

| File | Rows (committed) | Contents |
|---|---|---|
| `pbmc_baseline_links_full.csv` | 15,806 | Full retained LinkPeaks candidate set after `--candidate-filter positive_score` |
| `pbmc_baseline_links_with_distance.csv` | 15,806 | The above, plus peak/TSS coordinates, `distance_bp`, `distance_score` |
| `pbmc_link_features.csv` | 5,000 | Top `candidate_top_k` candidates with all scoring features. **This is the fixed evaluation universe.** |

Header of `pbmc_baseline_links_full.csv`:

```
peak,gene,seqnames,start,end,width,strand,score,zscore,pvalue
```

Header of `pbmc_link_features.csv` — 31 columns:

```
peak,gene,link_score,add,mul,mul_weigh,mul_strict,adj,
link_seqnames,link_start,link_end,link_width,link_strand,link_zscore,link_pvalue,
peak_chr,peak_start,peak_end,peak_mid,
gene_chr,tss,tx_id,tx_biotype,
distance_bp,distance_score,
peak_motif_score,peak_tf_score,peak_top_motifs,motif_score,tf_score,motif_names
```

### `results/<dataset>/scent_chr_sweep_<tag>/chr<N>/`

| File | Columns | Note |
|---|---|---|
| `scent_candidates_chr<N>.csv` | `gene,peak` | Cis-window candidate pairs generated by this repository's own code (`run_scent_chr_sweep.R` lines 325–365), not by SCENT. Across all 22 autosomes: **117,811 pairs, 9,891 genes, 46,936 peaks.** |
| `scent_links_chr<N>.csv` | `peak,gene,method,score,beta,se,z,p,boot_basic_p` | SCENT output |

`scent_chr_sweep_summary.csv`:

```
chr,status,runtime_minutes,n_rows,n_genes,n_peaks,n_positive,n_negative,
score_min,score_median,score_max
```

`status` takes values such as `skipped_existing` when `skip_existing: true` reused prior
output. Note that `runtime_minutes` is blank for skipped chromosomes, so the committed
summary does not record total compute time.

---

## 5. Final ranking outputs

One directory per score mode under `results/<dataset>/rankings/<mode>/`.

### `<dataset>_<mode>_ranked_links.csv` — the primary output

All 31 feature columns, plus 8:

| Column | Meaning |
|---|---|
| `model_score` | Final score $S_{pg}$ for this mode |
| `score_mode` | The `score_mode` value from `config/ablations.yaml` (note: `full_lambda_0_1` records `score_mode = full`, since the mode name and the score family differ) |
| `lambda_distance` | $\lambda$ used |
| `alpha_tf` | $\alpha$ used |
| `rank_model` | Rank by `model_score`, 1 = highest |
| `rank_link` | Rank by the LinkPeaks `link_score` on the same universe |
| `rank_diff_vs_linkpeaks` | `rank_link - rank_model`; positive = promoted by the model |
| `tier` | `High` / `Medium` / `Low`, from `tier_high_quantile` and `tier_medium_quantile` |

### Feature column meanings

| Column | Meaning |
|---|---|
| `peak` | Peak ID. **Format `chr1-1000-2000` (dashes).** See §9. |
| `gene` | Gene symbol |
| `link_score` | LinkPeaks score — the baseline ranking signal |
| `link_zscore`, `link_pvalue` | LinkPeaks statistics |
| `add`, `mul`, `mul_strict`, `adj` | Alternative coactivity formulations, retained for diagnostics; not used by any committed score mode |
| `mul_weigh` | **The coactivity score $A_{pg}$ used by all score modes.** Mean over cells of the product of positively-clipped RNA and ATAC z-scores |
| `peak_mid` | Peak midpoint, used for distance |
| `tss` | Representative TSS for the gene. One TSS per gene — see §9 |
| `tx_id`, `tx_biotype` | Transcript backing the chosen TSS |
| `distance_bp` | $\lvert$ `peak_mid` − `tss` $\rvert$ |
| `distance_score` | $D_{pg} = 1/(1 + (d/d_0)^2)$, with $d_0$ = `distance_d0` |
| `peak_motif_score` | Aggregate raw motif score for the peak |
| `peak_tf_score` | Aggregate motif × TF-expression score for the peak |
| `motif_score` | Scaled `peak_motif_score` |
| `tf_score` | **The TF/motif support score $T_p$ used by the score modes.** Scaled `peak_tf_score` |
| `peak_top_motifs`, `motif_names` | Up to `top_k_motif_names` motif identifiers and names |

`tf_score` is a **peak-level** quantity. It carries no gene-specific or cell-type-specific
information, which bounds what the TF term can contribute — see `docs/method_report.md` §7.

### Per-mode summaries

| File | Contents |
|---|---|
| `<dataset>_<mode>_summary_metrics.csv` | Long `metric,value` table: `n_links_ranked`, `n_unique_genes_ranked`, `n_full_linkpeaks_links_unrestricted`, `cor_linkpeaks_vs_model_score_same_universe`, top-50 median distances and distal fractions for model and baseline, ORA term counts, top-100 pair and gene overlap vs LinkPeaks |
| `<dataset>_<mode>_tier_summary.csv` | Per tier: `n_links`, `unique_genes`, `median_score`, `median_link_score`, `median_distance_bp` |
| `<dataset>_<mode>_top100_links.csv` | Top 100 by `model_score` |
| `<dataset>_<mode>_validation_topN_distance_diversity_summary.csv` | Top-N distance, distal fraction, gene diversity, link concentration |
| `<dataset>_<mode>_validation_component_correlations.csv` | Correlation of `model_score` with `link_score`, `mul_weigh`, `distance_score`, `tf_score` |
| `<dataset>_<mode>_validation_distance_matched_feature_contrast.csv` | Feature deltas, top decile vs rest, within distance bins |
| `<dataset>_<mode>_validation_distance_bin_rank_summary.csv` | Fuller within-bin rank summaries |
| `<dataset>_<mode>_validation_tf_motif_support_summary.csv` | TF/motif support across ranks, tiers, background |
| `<dataset>_<mode>_validation_promoted_demoted_inspection.csv` | Manual-review table for links promoted or demoted vs LinkPeaks |
| `<dataset>_<mode>_validation_manifest.csv` | Self-describing index of the six `validation_*` outputs. **Byte-identical across all modes** |

### Per-mode baseline artifacts — duplicated

These five files are re-emitted identically in every mode directory, because the LinkPeaks
baseline does not depend on the score mode:

```
<dataset>_<mode>_linkpeaks_baseline_same_universe.csv
<dataset>_<mode>_linkpeaks_baseline_distance_same_universe.csv
<dataset>_<mode>_linkpeaks_gene_rank_same_universe.csv
<dataset>_<mode>_linkpeaks_baseline_ora_GO_BP_same_universe.csv
<dataset>_<mode>_linkpeaks_baseline_ora_dotplot_same_universe.png
```

md5-verified identical across all 11 directories, totalling ~20.6 MB of redundancy.

### Plots

| File | Contents |
|---|---|
| `<dataset>_<mode>_distance_distribution.png` | Distance distribution of top-ranked links. The clearest single diagnostic for proximity collapse |
| `<dataset>_<mode>_linkpeaks_vs_model_scatter.png` | LinkPeaks score vs model score |
| `<dataset>_<mode>_ora_dotplot.png` | GO BP enrichment of top-N link genes |

---

## 6. SCENT validation outputs

`results/<dataset>/scent_validation/`

| File | Size | Contents |
|---|---|---|
| `scent_validation_scent_filter_summary.csv` | 191 B | One row: `scent_sweep_dir, support_rule, p_threshold, min_score, n_tested_rows, n_support_rows, n_tested_chromosomes` |
| `scent_validation_method_counts.csv` | 774 B | Per method: `n_links, unique_genes, unique_peaks, n_chromosomes, n_scent_supported, frac_scent_supported, score_{min,median,max}` |
| `scent_validation_topN_support_summary.csv` | 1.5 KB | **The headline table.** Per method × top-N: `n_scent_supported, frac_scent_supported, median_scent_rank_supported, median_distance_bp, distal_frac, promoter_frac_10kb, median_score` |
| `scent_validation_distance_matched_enrichment.csv` | 2.9 KB | Per distance bin × method: `n_high, n_rest, supported_high, supported_rest, frac_supported_{high,rest}, odds_ratio_high_vs_rest, median_distance_{high,rest}` |
| `scent_validation_supported_rank_summary.csv` | 570 B | Per method: rank position of SCENT-supported links |
| `scent_validation_pair_overlap_between_rankings.csv` | 4.5 KB | Pairwise top-K overlap between methods |
| `scent_validation_all_ranked_methods_combined.csv` | 6.1 MB | Long table, all methods × all candidates, with SCENT match columns. **Input to the min-distance summary** |
| `scent_validation_manifest.csv` | 779 B | Self-describing index |
| `top200_<method>_links_with_scent_support.csv` | ~45 KB × 6 | Top 200 per method, annotated |
| `top200_<method>_scent_supported_links.csv` | ~57 KB × 6 | Supported subset of the above |

### Interpreting `frac_scent_supported`

Two different quantities share this name:

- In `method_counts.csv` it is computed over the **whole** 4,976-pair universe. All six
  methods report **identically 1,259 / 4,976 = 0.2530**, because the universe is shared and
  only the ordering differs. This column carries no discriminative information.
- In `topN_support_summary.csv` it is computed over the **top N** of each method's ranking.
  This is the informative version.

### Columns of `scent_validation_all_ranked_methods_combined.csv`

```
peak,gene,method,score,original_score_mode,link_score,mul_weigh,distance_score,
distance_bp,tf_score,distance_modifier,rank_model,rank_link,rank,
scent_supported,scent_best_score,scent_best_rank,scent_best_peak,
scent_overlap_bp,scent_recip_overlap_query,scent_recip_overlap_scent,
scent_best_beta,scent_best_p,scent_best_boot_basic_p,scent_best_z
```

`scent_supported` is the boolean applied after the `pvalue_positive` rule. `scent_overlap_bp`
and the two reciprocal-overlap columns record how the ranked peak was matched to a SCENT
peak, controlled by `reciprocal_overlap: 0.5`.

### Plots

| File | Contents |
|---|---|
| `scent_validation_topK_supported_fraction.png` | SCENT-supported fraction vs rank depth, per method |
| `scent_validation_topK_median_distance.png` | Median distance vs rank depth. Read alongside the above — this is where `distance_only` is exposed |
| `scent_validation_distance_matched_enrichment_heatmap.png` | Odds ratios by distance bin × method |

---

## 7. Proximal-control outputs

`results/<dataset>/scent_validation_min_distance/`. Produced by
`summarize_scent_validation_min_distance.R`, consuming
`scent_validation_all_ranked_methods_combined.csv`, via `rule scent_validation_min_distance`
in `workflow/Snakefile` and `run_analysis.py run_scent_validation_min_distance`.

| File | Size | Contents |
|---|---|---|
| `scent_min_distance_method_counts.csv` | 1.3 KB | Per `min_distance_bp` × method: surviving `n_links`, `unique_genes`, `unique_peaks`, `n_scent_supported`, `frac_scent_supported_all`, `median_distance_bp` |
| `scent_min_distance_topN_support_summary.csv` | 4.2 KB | Per `min_distance_bp` × top-N × method, as the main top-N table |
| `scent_min_distance_delta_vs_linkpeaks.csv` | 4.9 KB | The above plus `linkpeaks_frac` and `delta_vs_linkpeaks`, sorted by support fraction |
| `scent_min_distance_distance_matched_enrichment.csv` | 8.2 KB | Distance-matched enrichment recomputed at each threshold |

Committed thresholds are 10,000 / 25,000 / 50,000 bp with top-N in {50, 100, 200, 500}.
The exact `--min-distances` and `--high-fraction` arguments used are not recorded in the
archive.

In `method_counts.csv`, the columns after `min_distance_bp` are **identical across all six
methods** at a given threshold — filtering removes the same pairs regardless of ranking.
Only `delta_vs_linkpeaks.csv` and `topN_support_summary.csv` distinguish methods.

No plot is produced. See `docs/results_report.md` §6 for the recommended figure and its path.

---

## 8. Commit / archive / exclude

### Safe to commit

| Path | Size |
|---|---|
| `config/`, `scripts/`, `workflow/`, `containers/`, `run_analysis.py` | < 250 KB |
| `renv.lock` | 90 KB |
| `docs/*.md` | ~120 KB |
| `results/pbmc/features/*.csv` | 7.4 MB |
| `results/pbmc/scent_validation/*.csv` except the combined table | ~330 KB |
| `results/pbmc/scent_validation/*.png` | 330 KB |
| `results/pbmc/scent_validation_min_distance/*.csv` | 36 KB |
| `results/pbmc/scent_chr_sweep_*/scent_chr_sweep_summary.csv` | ~3 KB |
| Per mode: `summary_metrics`, `tier_summary`, `top100_links`, `distance_distribution.png`, three `validation_*` summaries | ~140 KB × 11 |

Roughly **12–15 MB**, comfortable for GitHub.

### Should stay out of GitHub

| Path | Size | Reason |
|---|---|---|
| `data/` | large | Raw input; licensing unresolved |
| `resources/jaspar/JASPAR2022.sqlite` | large | Third-party database; redistribution terms unresolved |
| `renv/library/`, `renv/staging/` | very large | Restored inside the image |
| `tarballs/` | large | Nested duplicate of `results/` |
| `results/pbmc/rankings_backup/` | ~40 MB | Superseded 7-mode copy |
| `results/pbmc.before_restore_*/` | — | Restore scratch |
| `results/pbmc/rankings/*/pbmc_*_ranked_links.csv` | ~29 MB | Regenerable in minutes from the feature table |
| Duplicated baseline artifacts | ~20.6 MB | Ten redundant copies of five files |

### Candidates for Zenodo or Git LFS

These are now archived at
[10.5281/zenodo.22032568](https://doi.org/10.5281/zenodo.22032568) (CC BY 4.0), together with
`resources/jaspar/JASPAR2022.sqlite`. `scripts/download_inputs.sh` fetches the JASPAR file from
that deposit and verifies it against the tracked `.sha256`.


| Path | Size | Recommendation |
|---|---|---|
| All 11 `pbmc_*_ranked_links.csv` | ~29 MB | **Zenodo.** Fully regenerable |
| `scent_chr_sweep_*/chr*/scent_links_chr*.csv` (22) | ~15 MB | **Zenodo.** Expensive to regenerate — SCENT is the slow step |
| `scent_chr_sweep_*/chr*/scent_candidates_chr*.csv` (22) | ~2 MB | **Zenodo, and do not lose these.** The 117,811-pair de novo cis-window universe |
| `scent_links_all_chromosomes.csv` | 6.7 MB | **Zenodo** |
| `scent_validation_all_ranked_methods_combined.csv` | 6.1 MB | **Git LFS if available**, otherwise Zenodo. It has a live downstream consumer, so removing it from the clone breaks re-running the min-distance summary |
| `resources/jaspar/JASPAR2022.sqlite` | large | **LFS or runtime fetch**, subject to licensing |
| Per-mode `ora_dotplot.png`, `linkpeaks_vs_model_scatter.png` | ~7 MB | **Zenodo** |

---

## 9. Interface hazards

Four traps that will produce silently wrong results.

### 9.1 Peak ID format is inconsistent between subsystems

- LinkPeaks-derived files use **`chr11-60455251-60456168`** (all dashes)
- SCENT candidate files use **`chr1:229480033-229480809`** (colon, then dash)

`scripts/run_scent_chr_sweep.R` applies `normalize_peak()` internally, but files written by
the two halves of the pipeline are **not** in the same format on disk. Joining
`pbmc_baseline_links_full.csv` to `scent_candidates_chr*.csv` on the raw `peak` string
returns **zero** matches. After normalisation the true overlap is 5,768 pairs.

This hazard is still live. Always split on `[:-]` before comparing peak identifiers across
subsystems.

### 9.2 The evaluated universe is 4,976, not 5,000

`restrict_to_scent_chrs: true` drops candidates on contigs the SCENT sweep did not cover.
Every SCENT-validation number is over 4,976 pairs, 1,375 genes, 4,087 peaks. The 5,000 /
1,390 figures describe the feature table before restriction. Do not mix them.

### 9.3 Candidate window and validation window differ

LinkPeaks candidates use `link_distance: 500000`. The SCENT sweep uses `link_distance: 100000`.
SCENT support is therefore **structurally unavailable** for any candidate beyond 100 kb.
Distance bins at and above 100 kb cannot be interpreted as negative evidence — they are
untested. This drives the degenerate rows discussed in `docs/results_report.md` §5.

### 9.4 SCENT was run without cell-type stratification

`scoring_celltype: ""` means a synthetic `all_cells` label. Nothing in the committed results
is cell-type-specific, on either the ranking or the validation side.

---

## 10. Reproducing each stage

Via the wrapper (defaults `--image multiome-reranking-benchmark:v0.1.0`,
`--snakefile workflow/Snakefile`, `--configfile config/default.yaml`, `--cores 4`):

```bash
python3 run_analysis.py list_score_modes
python3 run_analysis.py build_linkpeaks_features
python3 run_analysis.py run_score_mode --mode full_lambda_0_1
python3 run_analysis.py run_all_score_modes
python3 run_analysis.py run_scent_sweep
python3 run_analysis.py run_scent_validation
python3 run_analysis.py run_scent_pipeline
python3 run_analysis.py run_all_score_modes --dry-run
python3 run_analysis.py unlock
```

Directly:

```bash
snakemake --snakefile workflow/Snakefile --configfile config/default.yaml --cores 4 all
snakemake --snakefile workflow/Snakefile --configfile config/default.yaml --cores 4 all_with_scent
```

The min-distance controls currently have no target and must be run by hand:

```bash
Rscript scripts/summarize_scent_validation_min_distance.R \
  --input results/pbmc/scent_validation/scent_validation_all_ranked_methods_combined.csv \
  --output-dir results/pbmc/scent_validation_min_distance \
  --min-distances 10000,25000,50000 \
  --top-n-values 50,100,200,500
```

`--high-fraction` was not recorded and is left at its script default.
