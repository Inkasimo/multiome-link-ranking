# Method report — PBMC LinkPeaks-candidate peak–gene reranking benchmark

This document describes what the current repository implements. It is a description of a
**reranking benchmark over a fixed candidate universe**, not of a standalone peak–gene
linking method. Future work is described separately in `docs/future_standalone_v0.md`.

All formulas below are transcribed from the implementation, with file and line references.
Where the implementation differs from the conceptual description that circulated during
development, the implementation is authoritative and the difference is noted.

---

## 1. Question addressed

The benchmark addresses a narrow, deliberately answerable question:

> Given a fixed set of candidate peak–gene pairs produced by `Signac::LinkPeaks()`, does
> reordering those same pairs using interpretable RNA–ATAC coactivity, genomic distance, and
> TF/motif support produce a ranking with more external support than the LinkPeaks ordering
> itself?

Three properties of this framing are load-bearing:

1. **The candidate universe is held fixed.** Every score mode ranks the identical set of
   pairs. No mode can gain by admitting or excluding candidates.
2. **The comparison is ordinal.** No calibrated probability or effect size is produced.
3. **The baseline is also the candidate generator.** LinkPeaks both defines the universe and
   supplies the reference ordering. This is a structural constraint on interpretation, not a
   bug — see §12.

The question deliberately excluded from scope: whether the scoring components can generate a
candidate universe of their own. That requires de novo candidate generation and is not
implemented here.

---

## 2. Input data assumptions

| Assumption | Basis |
|---|---|
| Paired single-cell RNA + ATAC from the same cells (10x multiome) | `Read10X_h5()` on a filtered feature-barcode matrix carrying both `Gene Expression` and `Peaks` assays |
| Dataset: 10x `pbmc_unsorted_10k`, Cell Ranger ARC 2.0.0, reference `refdata-cellranger-arc-GRCh38-2020-A-2.0.0` | Recovered from the `atac_fragments.tsv.gz` header (`id=`, `pipeline_version=`, `reference_path=`, `reference_version=2020-A`) |
| ~12,012 cells, 111,857 called ATAC peaks, one healthy female donor | 10x dataset metrics; 50,918 peaks (45.5%) survive into the committed results |
| Human, hg38 | `genome(annotations) <- "hg38"` (`run_linkpeaks_reranker.R:377`), `BSgenome.Hsapiens.UCSC.hg38` (`:446`) |
| Ensembl v86 annotation | `GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)` (`:375`) |
| Tabix-indexed fragments | `.tbi` declared as an explicit Snakemake input |
| Vertebrate JASPAR2022 CORE motifs, human | `getMatrixSet(JASPAR2022, list(collection = "CORE", tax_group = "vertebrates", species = 9606))` (`:616`) |
| Single sample, no batch structure | No integration or batch-correction step exists in the pipeline |
| No cell-type annotation | Clustering is performed at `cluster_resolution: 0.5` on the WNN graph; clusters are never labelled, and **no score in this repository is cell-type-stratified** |

**Dataset provenance was not recorded in the repository and had to be recovered forensically**
from the `cellranger-arc` comment block at the head of `data/atac_fragments.tsv.gz`:

```
# id=pbmc_unsorted_10k
# description=PBMC from a healthy donor - no cell sorting (10k)
# pipeline_name=cellranger-arc
# pipeline_version=cellranger-arc-2.0.0
# reference_path=.../refdata-cellranger-arc-GRCh38-2020-A-2.0.0
# reference_fasta_hash=b6f131840f9f337e7b858c3d1e89d7ce0321b243
# reference_gtf_hash=3b4c36ca3bade222a5b53394e8c07a18db7ebb11
# reference_version=2020-A
```

This is 10x Genomics' publicly available `pbmc_unsorted_10k` dataset (published 2021-05-03,
CC BY 4.0). It should be recorded in `config/default.yaml` so it does not have to be recovered
again — see `TODO.md` §8.

Note this is **not** the `pbmc_granulocyte_sorted_10k` dataset used in the Signac multiome
vignette, despite the identical local filenames. The fragment-file sizes differ
(2,917,757,251 B here versus 2,051,027,831 B for the sorted dataset), which is an independent
confirmation.

### Preprocessing

**No cell-level quality control is performed.** `CreateSeuratObject()` is called at
`run_linkpeaks_reranker.R:380` without `min.cells` or `min.features`, and there is no `subset()`
step anywhere in the script. All barcodes present in `filtered_feature_bc_matrix.h5` — roughly
12,012 for this dataset — enter the analysis. Comparable Signac multiome workflows filter on
`nCount_ATAC`, `nCount_RNA`, `nucleosome_signal` and `TSS.enrichment` before proceeding, which
typically removes several percent of barcodes. Low-quality nuclei therefore contribute to every
per-cell z-score and to \(A_{pg}\).

1. RNA: normalisation, variable-feature selection, scaling, PCA to `pca_dims: 30`.
2. ATAC: TF-IDF, SVD, LSI components `lsi_dims_start: 2` through `lsi_dims_end: 30`. Dropping
   component 1 is standard, as it correlates with sequencing depth.
3. WNN integration of the two modalities, then Louvain clustering.
4. `RegionStats()` against hg38, required by `LinkPeaks()` for GC and sequence covariates.

Clustering output is not consumed by any scoring step. It exists because `LinkPeaks()` is run
on a fully processed object, not because the score uses cluster identity.

---

## 3. Candidate universe

Two nested sets:

| Set | Size | Definition |
|---|---|---|
| Retained LinkPeaks candidates | **15,806** | `LinkPeaks()` output at `link_distance: 500000`, filtered by `--candidate-filter positive_score` |
| Feature-table universe | **5,000** | Top `candidate_top_k: 5000` of the above |
| SCENT-evaluable subset | **4,976** | The 5,000 restricted to chromosomes the SCENT sweep covered (`restrict_to_scent_chrs: true`) |

The 5,000-pair table spans **1,390 genes**. The 4,976-pair evaluated subset spans **1,375
genes and 4,087 peaks**.

The distinction between 5,000 and 4,976 matters and is easy to lose. Every number in the
SCENT validation is over 4,976. Every number in the per-mode ranking diagnostics is over 5,000.

`--candidate-filter` is not exposed in `config/default.yaml`, so the filter defining the
universe is currently an invisible script default.

---

## 4. LinkPeaks baseline

`Signac::LinkPeaks()` is used in two roles simultaneously:

1. **Candidate generator** — produces the pairs.
2. **Baseline ranking** — its `score` column becomes `link_score`, and the `linkpeaks` score
   mode ranks by it directly (`evaluate_rankings.R:489`).

LinkPeaks correlates peak accessibility with gene expression across cells and calibrates
against a background matched on GC content, accessibility, and peak width, yielding a z-score
and p-value. The `positive_score` filter retains pairs with positive association.

Because the baseline is also the candidate generator, "beats LinkPeaks" means "reorders
LinkPeaks' own output more usefully than LinkPeaks does", not "finds links LinkPeaks missed".
Recall is bounded by LinkPeaks at every mode.

---

## 5. Feature table construction

Built once by `scripts/run_linkpeaks_reranker.R`; all eleven committed score modes read from it. This is
the central engineering decision — the heavy Seurat/Signac/motif work is not repeated per mode.

Stages:

1. Build the Seurat object, preprocess both modalities, integrate, cluster.
2. Run `LinkPeaks()`; apply the candidate filter; take the top `candidate_top_k`.
3. Compute per-cell z-scores for each candidate gene and peak.
4. Compute coactivity variants (§6).
5. Attach peak midpoints and gene TSS; compute distance and distance score (§7).
6. Match JASPAR motifs to peaks; compute peak-level motif and TF scores (§8).
7. Write the three CSVs in `results/<dataset>/features/`.

Five coactivity variants are computed and stored (`add`, `mul`, `mul_weigh`, `mul_strict`,
`adj`). Only `mul_weigh` is used by any committed score mode; the others are retained for
diagnostics.

---

## 6. Coactivity score

For candidate pair \((p, g)\) over cells \(c \in C\), let \(z^{\mathrm{RNA}}_{gc}\) and
\(z^{\mathrm{ATAC}}_{pc}\) be the z-scored RNA expression of \(g\) and accessibility of \(p\).

The coactivity score used throughout is `mul_weigh` (`run_linkpeaks_reranker.R:554`):

$$
A_{pg} \;=\; \frac{1}{|C|} \sum_{c \in C} \max\!\left(z^{\mathrm{RNA}}_{gc},\, 0\right) \cdot \max\!\left(z^{\mathrm{ATAC}}_{pc},\, 0\right)
$$

Only cells above the mean in **both** modalities contribute. The positive clipping makes this
a soft AND: a pair is rewarded when accessibility and expression are jointly elevated in the
same cells, and cells where either is below its mean contribute exactly zero rather than a
negative term.

The other stored variants, for reference:

$$
\mathrm{add}_{pg} = \frac{1}{|C|}\sum_c \left[\max(z^{\mathrm{RNA}}_{gc},0) + \max(z^{\mathrm{ATAC}}_{pc},0)\right]
\qquad
\mathrm{mul}_{pg} = \frac{1}{|C|}\sum_c z^{\mathrm{RNA}}_{gc} z^{\mathrm{ATAC}}_{pc}
$$

$$
\mathrm{mul\_strict}_{pg} = \frac{1}{|C|}\sum_c \mathbb{1}\!\left[z^{\mathrm{RNA}}_{gc} > 1\right]\mathbb{1}\!\left[z^{\mathrm{ATAC}}_{pc} > 1\right]
\qquad
\mathrm{adj}_{pg} = \frac{\mathrm{mul\_strict}_{pg}}{\bar{r}_g \bar{a}_p + 10^{-6}}
$$

where \(\bar{r}_g = \frac{1}{|C|}\sum_c \mathbb{1}[z^{\mathrm{RNA}}_{gc} > 1]\) and
\(\bar{a}_p\) is its ATAC analogue.

Two properties of \(A_{pg}\) constrain interpretation. It is computed over **all cells
pooled**, so a link active in one small population is diluted. And because both inputs are
z-scored per feature, \(A_{pg}\) is scale-free but its absolute magnitude is not comparable
across datasets.

---

**\(A_{pg}\) is associated with marginal detection rates.** Recovering the marginal activity
product as `mul_strict / adj` gives Spearman(`mul_weigh`, marginals) = **+0.682**. `mul_weigh` has
no natural null: under independence its expectation is positive and depends on each feature's
detection rate. LinkPeaks controls for this via a background matched on accessibility and GC;
`mul_weigh` does not. **This is not conditioned on anywhere in the benchmark.**

The feature table contains an adjusted variant, `adj` = `mul_strict` / (marginal activity product)
— a lift ratio, null at 1. It is stored but **not used by any score mode**, and was inspected in
early prototype work without being adopted. Re-examination shows the ratio normalisation
overcorrects: Spearman(`adj`, marginals) = **−0.867**, Spearman(`adj`, `mul_strict`) = **−0.735**.
Dividing by a denominator 0.968-collinear with the numerator yields a ratio governed by the
denominator, so `adj` preferentially ranks rare gene × rare peak pairs. Recorded as a failed
correction, not a robustness result.

## 7. Distance score

The gene TSS is taken per candidate pair as the transcript TSS **minimising** the distance to
the peak midpoint (`run_linkpeaks_reranker.R:305–316`). So

$$
d_{pg} \;=\; \min_{t \in \mathcal{T}(g)} \left| m_p - \mathrm{tss}_t \right|
$$

where \(m_p\) is the peak midpoint and \(\mathcal{T}(g)\) the transcripts of \(g\). The
transcript backing the selection is recorded in `tx_id` / `tx_biotype`.

The distance score is a squared-Lorentzian decay (`:319–322`):

$$
D_{pg} =
\begin{cases}
\frac{1}{1 + (d_{pg}/d_0)^2}, & d_{pg}\ \mathrm{finite} \\[2ex]
0, & \mathrm{otherwise}
\end{cases}
$$

where \(d_0 = 50{,}000\), from config field `distance_d0`.

so \(D \to 1\) at the TSS, \(D = 0.5\) at \(d = d_0\), and \(D \approx 0.01\) at 500 kb. The
tail is heavy — distal links are downweighted, never excluded.

**A TSS-convention inconsistency exists between subsystems.** The reranker takes the
closest transcript TSS per pair, as above. `run_scent_chr_sweep.R:138` uses "one
representative TSS per gene". Distances computed by the two halves of the pipeline are
therefore not guaranteed identical for the same pair. This does not affect the internal
consistency of the reranking comparison, since all modes share one feature table, but it does
affect distance-stratified comparisons against SCENT.

---

## 8. TF/motif support score

Computed per **peak**, then attached to every pair containing that peak
(`run_linkpeaks_reranker.R:616–695`).

1. Retrieve JASPAR2022 CORE vertebrate PFMs for species 9606.
2. Retain motifs whose associated TF is expressed in at least
   `tf_expressed_frac: 0.10` of cells.
3. Run `motifmatchr` against hg38 to obtain a continuous peak × motif score matrix
   \(M \in \mathbb{R}^{P \times K}\), column-wise min–max rescaled (`:647`).
4. For each motif \(k\), form a TF-expression weight from the **global mean expression** of
   its associated TF genes, taking the maximum over TFs mapped to that motif, then min–max
   rescale across motifs (`:657–664`):

$$
w_k \;=\; \mathrm{rescale}_{01}\!\left( \max_{j \in \mathcal{G}(k)} \; \overline{x}^{\mathrm{RNA}}_{j} \right),
\qquad \overline{x}^{\mathrm{RNA}}_{j} = \frac{1}{|C|}\sum_{c} x^{\mathrm{RNA}}_{jc}
$$

5. Combine and rescale (`:666–667`):

$$
T_p \;=\; \mathrm{rescale}_{01}\!\left( \sum_{k=1}^{K} M_{pk} \, w_k \right) \;\in\; [0, 1]
$$

with \(\mathrm{rescale}_{01}(x) = (x - \min x)/(\max x - \min x)\), returning all zeros for a
degenerate range (`:145–155`).

Three consequences follow directly and bound what this term can do:

- **\(T_p\) has no gene dependence.** All pairs sharing a peak receive the same \(T_p\). The
  term cannot express "this TF regulates *this* gene"; it expresses "this peak looks like a
  regulatory element".
- **\(T_p\) has no cell-type dependence.** \(w_k\) uses global mean expression.
- **\(T_p\) is min–max scaled within the dataset**, so it is not portable across datasets and
  is sensitive to outliers at either end.

`motif_score` (= `peak_motif_score`) is the same construction without the TF-expression
weighting, and is stored but unused by the committed modes.

---

## 9. Score family

All modes share one multiplicative template:

```text
S_pg = A_pg * f_lambda(D_pg) * g_alpha(T_p)

g_alpha(T_p) = 1 + alpha * T_p
alpha = alpha_tf
```

Multiplicative rather than additive composition is intentional: weak evidence in any one
component should reduce confidence rather than be offset by strength elsewhere. Because
\(g_\alpha \geq 1\) and \(f_\lambda > 0\), neither modifier can zero out a link — they
reweight rather than veto.

### 9.1 Original distance prior

Implemented for `coactivity_distance`, `full`, and `full_linkpeaks_anchored`
(`evaluate_rankings.R:502, 514`):

$$
f^{\mathrm{orig}}_{\lambda}(D) \;=\; (1 - \lambda) + \lambda D
$$

Range \([1-\lambda,\ 1]\). At \(\lambda = 0.1\): \([0.90,\ 1.00]\). The modifier is
**always ≤ 1**, so this form is a pure penalty on distance — no link is boosted relative to
an unmodified coactivity ranking, and the most TSS-proximal link is merely left untouched.

### 9.2 Modified distance prior

Implemented for `distance_mod_only` and `full_moddist` (`evaluate_rankings.R:79–84`):

$$
f^{\mathrm{mod}}_{\lambda}(D) \;=\; 1 + \lambda\left(D - \tfrac{1}{2}\right)
$$

Range \([1 - \lambda/2,\ 1 + \lambda/2]\). At \(\lambda = 0.1\): \([0.95,\ 1.05]\). This is
symmetric about \(D = 0.5\), i.e. about \(d_{pg} = d_0 = 50\) kb: links closer than 50 kb are
mildly boosted, links further are mildly penalised, and a link at exactly 50 kb is unchanged.

Conceptually this is the cleaner object — it separates "distance is informative" from
"distance is a penalty", and the neutral point is an interpretable quantity rather than an
artifact of the parameterisation.

The two forms are related by an affine map:

$$
f^{\mathrm{mod}}_{\lambda}(D) \;=\; f^{\mathrm{orig}}_{\lambda}(D) + \tfrac{\lambda}{2}
$$

The offset \(\lambda/2\) is **constant across pairs**. Multiplying \(A_{pg} g_\alpha(T_p)\)
by a constant-shifted modifier does not preserve rank order in general — the shift interacts
with the magnitude of \(A_{pg} g_\alpha(T_p)\) — but the induced reordering is small, which is
consistent with the observed near-identity of the two families (§ `docs/results_report.md` 7).

---

## 10. Score modes and ablations

Twelve `score_mode` branches are implemented in `evaluate_rankings.R`; eleven are configured in
`config/ablations.yaml`, and all eleven have committed results. The twelfth,
`full_linkpeaks_anchored` (`evaluate_rankings.R:527`), has no configuration entry and was never
run.

| Mode | `score_mode` | \(\lambda\) | \(\alpha\) | \(S_{pg}\) | Role |
|---|---|---|---|---|---|
| `linkpeaks` | `linkpeaks` | — | — | `link_score` | Baseline |
| `coactivity` | `coactivity` | — | — | \(A_{pg}\) | Isolate RNA–ATAC evidence |
| `distance_only` | `distance_only` | — | — | \(D_{pg}\) | **Negative control** |
| `distance_mod_only_lambda_0_1` | `distance_mod_only` | 0.10 | 0.00 | \(f^{\mathrm{mod}}_{0.1}(D_{pg})\) | Control on the modified form |
| `coactivity_distance` | `coactivity_distance` | 0.30 | 0.00 | \(A \cdot f^{\mathrm{orig}}_{0.3}\) | Add distance |
| `coactivity_tf` | `coactivity_tf` | 0.00 | 0.50 | \(A \cdot (1 + 0.5\,T)\) | Add TF/motif |
| `full` | `full` | 0.30 | 0.50 | \(A \cdot f^{\mathrm{orig}}_{0.3} \cdot (1 + 0.5\,T)\) | Full, strong distance |
| `full_lambda_0_1` | `full` | 0.10 | 0.50 | \(A \cdot f^{\mathrm{orig}}_{0.1} \cdot (1 + 0.5\,T)\) | Full, mild distance |
| `full_lambda_0_2` | `full` | 0.20 | 0.50 | \(A \cdot f^{\mathrm{orig}}_{0.2} \cdot (1 + 0.5\,T)\) | Full, intermediate |
| `full_moddist_lambda_0_1` | `full_moddist` | 0.10 | 0.50 | \(A \cdot f^{\mathrm{mod}}_{0.1} \cdot (1 + 0.5\,T)\) | Modified distance, mild |
| `full_moddist_lambda_0_2` | `full_moddist` | 0.20 | 0.50 | \(A \cdot f^{\mathrm{mod}}_{0.2} \cdot (1 + 0.5\,T)\) | Modified distance, intermediate |

Note that the mode **name** and the `score_mode` field differ: `full_lambda_0_1` records
`score_mode = full` in its output. Read the `lambda_distance` and `alpha_tf` columns to
disambiguate.

`distance_only` is the important one. It is included not as a competitor but as a detector:
if a proposed score cannot outperform pure TSS proximity, the score is measuring proximity.
Its behaviour is the primary interpretive constraint on the whole benchmark.

### Ablation logic

- `coactivity` vs `linkpeaks` — does raw coactivity carry independent signal?
- `coactivity_distance` vs `coactivity` — does the distance prior add or merely reshuffle?
- `coactivity_tf` vs `coactivity` — does peak-level TF/motif support add?
- `full` vs its parts — do the components combine constructively?
- `full_lambda_0_1`, `full_lambda_0_2`, `full` (λ = 0.3) — sensitivity to \(\lambda\).
- `full_moddist_*` vs `full_*` — does the reparameterisation change anything empirically?
- `distance_only`, `distance_mod_only_*` — controls.

Seven of the eleven modes were carried into the SCENT comparison
(`config/scent_validation.yaml`): `linkpeaks`, `coactivity`, `coactivity_tf`,
`full_lambda_0_1`, `full_moddist_lambda_0_1`, `full`, `distance_only`. Four —
`coactivity_distance`, `full_lambda_0_2`, `full_moddist_lambda_0_2` and
`distance_mod_only_lambda_0_1` — have no external comparison.

`full` (λ = 0.3) is included as an aggressive distance-prior sensitivity setting, not as a
candidate primary model. `full_lambda_0_1` (λ = 0.1) remains the conservative primary setting.

SCENT is a comparator, not ground truth: it consumes the same RNA and ATAC matrices, is
correlational, and was run with a 100 kb candidate window while these candidates extend to
500 kb.

---

## 11. SCENT validation

### 11.1 Rationale

Gene-level over-representation analysis proved insufficient. ORA evaluates the *gene set*
implied by top-ranked links and is blind to which peak was paired with which gene — a
ranking could be scrambled at the peak level and yield an identical ORA result. A link-level
comparator was required.

SCENT (`immunogenomics/SCENT` v1.0.1, commit `e80b5ba6b445f972c7fe28fb41e24ef4f5b2e373`)
fits a Poisson model of gene counts on peak accessibility with bootstrap inference, per
peak–gene pair. It is used here as an **independent correlational comparator**, not as ground
truth. It shares both input matrices with the reranker, so agreement is not independent
validation in the strong sense.

### 11.2 Candidate generation for the sweep

`scripts/run_scent_chr_sweep.R` builds its own candidate set — this code is part of this
repository, not part of SCENT:

1. `expr_frac_gene <- Matrix::rowMeans(rna_counts > 0)`,
   `expr_frac_peak <- Matrix::rowMeans(atac_counts > 0)` (`:325–326`)
2. Retain genes and peaks with fraction \(\geq\) `min_pair_frac: 0.02` (`:328–329`)
3. Same-chromosome cartesian join of gene TSS × peaks (`:353–358`)
4. Retain pairs with \(\lvert m_p - \mathrm{tss}_g \rvert \leq\) `link_distance: 100000` (`:364–365`)

Across all 22 autosomes this produced **117,811 candidate pairs over 9,891 genes and 46,936
peaks** — a cis-window universe independent of LinkPeaks. Overlap with the 15,806 LinkPeaks
candidates, after peak-ID normalisation, is 5,768 pairs (36.5% of the LinkPeaks set; 42.2%
restricting to shared genes).

Execution used `max_cells: 1000`, `scent_regr: poisson`, `scent_cores: 4`, and
`scoring_celltype: ""` — meaning a synthetic `all_cells` label, so **the SCENT results are
not cell-type-stratified**.

Of the 117,811 candidates, **52,482 rows** were tested and returned by SCENT across 22
autosomes.

### 11.3 Support rule

`scent_support_rule: pvalue_positive` with `scent_p_threshold: 0.05`. A SCENT row is
supporting if

$$
\beta > 0 \quad\text{and}\quad p_{\mathrm{boot}} \leq 0.05
$$

Applied to the 52,482 tested rows, this retains **4,758 support rows**.

### 11.4 Matching ranked links to SCENT rows

Peak intervals differ between subsystems, so matching is by reciprocal overlap. A ranked pair
\((p, g)\) matches a SCENT pair \((p', g)\) when the gene is identical and

$$
\frac{\lvert p \cap p' \rvert}{\lvert p \rvert} \geq 0.5
\quad\text{and}\quad
\frac{\lvert p \cap p' \rvert}{\lvert p' \rvert} \geq 0.5
$$

(`reciprocal_overlap: 0.5`). The best-scoring match is retained, and the overlap statistics
are written to `scent_overlap_bp`, `scent_recip_overlap_query`, `scent_recip_overlap_scent`.

With `restrict_to_scent_chrs: true`, ranked methods are cut to SCENT-covered chromosomes,
giving the **4,976-pair** evaluated universe.

### 11.5 Metrics

- **Top-N supported fraction** at N ∈ {50, 100, 200}: of a method's top N, the fraction
  matching a supporting SCENT row. The headline metric.
- **Rank of supported links**: median and mean `rank_model` over supported pairs. Lower is
  better; unlike the top-N fraction this uses the full ranking.
- **Median distance and distal fraction** at each N, plus `promoter_frac_10kb`. Reported
  *alongside* the support fraction, because a method can raise support by becoming more
  proximal.
- **Pairwise top-K overlap** between methods, at `top_k_compare: 200`.

---

## 12. Distance-matched validation

The core confound: SCENT support is itself strongly enriched near promoters, and the distance
score is a monotone function of proximity. A method can therefore improve its support
fraction purely by ranking proximal links higher, with no additional biological insight.

The mitigation is stratification. Candidates are binned by \(d_{pg}\) into
`0_10kb`, `10_50kb`, `50_200kb`, `200_500kb`, `gt500kb`. Within each bin, the method's top
decile is compared against the remainder:

$$
\mathrm{OR}_{\text{bin}} \;=\;
\frac{\mathrm{supported\_high} \,/\, (n_{\mathrm{high}} - \mathrm{supported\_high})}
     {\mathrm{supported\_rest} \,/\, (n_{\mathrm{rest}} - \mathrm{supported\_rest})}
$$

An odds ratio above 1 within a bin indicates the score orders links usefully *at fixed
distance*, which proximity alone cannot explain.

Two limits on this analysis, both structural:

1. **Bin occupancy is fixed by the candidate set, not by the method.** `n_high` and `n_rest`
   are identical across methods within a bin (150 / 1,349 in `0_10kb`; 111 / 999 in
   `10_50kb`). Only the supported counts vary.
2. **Bins beyond 100 kb are untested, not negative.** The SCENT sweep used a 100 kb window
   while candidates extend to 500 kb. In `200_500kb` and `gt500kb`, `supported_high` and
   `supported_rest` are both **0** for every method — no SCENT test exists there. The odds
   ratios reported for those bins (8.906 and 0.333) are continuity-correction artifacts on
   empty cells and carry no information. The `50_200kb` bin is partially affected: only its
   50–100 kb portion is testable.

Interpretable bins are therefore `0_10kb`, `10_50kb`, and — with the above caveat —
`50_200kb`.

---

## 13. Proximal and min-distance controls

A second, blunter control, implemented in
`scripts/summarize_scent_validation_min_distance.R`. Rather than stratifying by distance, it
**removes** all links below a threshold and recomputes the comparison on what remains:

It is wired into the workflow as `rule scent_validation_min_distance` (`workflow/Snakefile`),
configured by `config/scent_validation_min_distance.yaml`, and invoked as
`python3 run_analysis.py run_scent_validation_min_distance`. It is post-processing of the
validation output and does not re-run SCENT.

```text
C_delta = {(p, g): d_pg > delta}
delta in {10 kb, 25 kb, 50 kb}
```

Top-N support fractions are recomputed within \(\mathcal{C}_\delta\) at
N ∈ {50, 100, 200, 500}, and `delta_vs_linkpeaks` records each method's advantage over the
baseline at that threshold.

The logic: if a method's apparent advantage is carried by promoter-proximal links, removing
them should erase it.

Filtering is method-independent — \(\mathcal{C}_\delta\) is the same set for all methods —
so `scent_min_distance_method_counts.csv` shows identical surviving counts across methods at
each threshold (3,477 at 10 kb; 2,871 at 25 kb; 2,367 at 50 kb).

**This control has a ceiling that must be stated.** As \(\delta\) grows, the surviving median
distance grows: 105.8 kb at \(\delta\) = 10 kb, 150.0 kb at 25 kb, 196.2 kb at 50 kb. But
SCENT only tested pairs within 100 kb. At \(\delta\) = 50 kb the surviving candidates are
mostly beyond SCENT's tested range, so support can only arise from the narrow 50–100 kb band.
The control becomes progressively less able to discriminate exactly as it becomes most
relevant. The \(\delta\) = 50 kb result should be read as weak evidence in either direction,
not as a clean negative.

The overall SCENT support rate falls sharply with \(\delta\) — 0.177 at 10 kb, 0.134 at
25 kb, 0.083 at 50 kb — which is consistent both with genuine biology (proximal links are
more often real) and with the window artifact. These two explanations cannot be separated
using SCENT alone.

---

## 14. Limitations

Ordered by how much they constrain the conclusions.

### Structural

1. **LinkPeaks defines the candidate universe.** No mode can exceed LinkPeaks' recall. The
   benchmark measures reordering quality only.
2. **The validator is not ground truth.** SCENT consumes the same two matrices as the
   reranker and is itself correlational, promoter-biased, and window-limited. Concordance
   between two correlational methods on shared input is weaker evidence than it appears.
3. **The validation window is narrower than the candidate window** — 100 kb vs 500 kb. The
   distal regime, which is the field's actual open problem, is unmeasured here.
4. **`distance_only` is not cleanly beaten.** It wins at top-50 on the full universe and
   remains the strongest method at every proximal-removal threshold. See
   `docs/results_report.md` §5, §6.

### Methodological

5. **No cell-type stratification anywhere.** \(A_{pg}\) pools all cells; \(w_k\) uses global
   mean expression; SCENT ran with a synthetic `all_cells` label. Links active in a single
   population are systematically diluted.
6. **\(T_p\) is peak-level only.** It cannot represent TF-to-target specificity, which is
   what "TF support for this link" would require.
7. **Coactivity is not conditioned on marginal activity.** `mul_weigh` correlates with the
   marginal detection-rate product at +0.682. Whether its contribution survives conditioning on
   marginal activity **is not tested in this release**. The natural check — activity-matched
   enrichment mirroring §12 — is deferred to v0.2.
8. **Min–max rescaling of \(T_p\), \(w_k\), and the motif matrix** makes scores
   dataset-dependent and outlier-sensitive. No score in this repository is portable across
   datasets.
9. **TSS conventions differ between subsystems** (§7): closest transcript TSS in the
   reranker, one representative TSS in the SCENT sweep.
10. **\(\lambda\) and \(\alpha\) are set by hand, not fitted.** Deliberate, for
   interpretability, but it means the ablation grid is a sensitivity sweep and not an
   optimisation. No held-out selection was performed.
11. **\(f^{\mathrm{orig}}_\lambda \leq 1\) always** (§9.1), so the original distance prior can
    only demote. This asymmetry was unintended and motivated the modified form.

### Scale and generality

12. **One dataset, one tissue, one sample.** No replicate, no second tissue, no cross-dataset
    check.
13. **5,000 pairs over 1,390 genes** — a small fraction of the cis-regulatory space, and
    already the top of a LinkPeaks ranking, so not a random sample of candidates.
14. **`max_cells: 1000` for SCENT** but all cells for the reranker. The two sides of the
    comparison were computed on different cell sets.
15. **No orthogonal validation.** No CRISPRi perturbation data, no fine-mapped eQTLs, no
    chromatin-contact data. All available evidence is derived from the same two matrices.

### Data provenance

16. **No cell-level QC.** All ~12,012 barcodes are used, with no filtering on counts,
    nucleosome signal or TSS enrichment. Low-quality nuclei contribute to every coactivity
    estimate.
17. **Annotation versions are mismatched between quantification and TSS assignment.** The count
    matrix comes from the Cell Ranger ARC `GRCh38-2020-A` reference (a GENCODE-based build; the
    exact GENCODE/Ensembl correspondence should be confirmed from 10x's reference build notes
    before citing). TSS coordinates driving `distance_bp` and `distance_score` come from
    `EnsDb.Hsapiens.v86`, i.e. Ensembl 86 (2016). Genes present in the matrix but absent from
    Ensembl 86 are silently dropped, and TSS positions may differ for genes present in both.
    Since distance is the central confound in this benchmark, this is worth stating explicitly.
    The Signac multiome vignette pairs the same two annotation sources, so this is a widespread
    convention rather than an unusual error — but it is a real source of coordinate drift.

### Reproducibility

18. **Dataset provenance is resolved.** Dataset provenance, download URLs, file sizes and
    SHA256 checksums are now recorded in `config/default.yaml`.
19. **`--candidate-filter` not exposed** in configuration.
20. ~~The min-distance controls have no Snakemake rule, and the exact arguments used for the
    committed outputs are not recorded.~~ **Resolved.** Implemented as
    `rule scent_validation_min_distance` (`workflow/Snakefile`), configured by
    `config/scent_validation_min_distance.yaml` (`min_distances: 10000,25000,50000`,
    `top_n_values: 50,100,200,500`, `high_fraction: 0.10`), exposed as
    `python3 run_analysis.py run_scent_validation_min_distance`.
21. **The JASPAR sqlite pinning is undocumented** outside the script body, and the file is
    not in the container image.
