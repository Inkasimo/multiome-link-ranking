# Positioning against similar tools

Honest positioning of this repository, and of the future standalone direction, against existing
peak–gene and cis-regulatory inference methods.

The purpose of this document is defensive, not promotional. It records which claims are
unavailable and which narrower framings survive scrutiny. Two things should be stated at the top:

- **The current repository is a benchmark, not a method.** It reranks a `Signac::LinkPeaks()`
  candidate universe and evaluates the result. Against every method below, its correct
  positioning is "downstream evaluation of one of these, using another as a comparator", not
  "alternative to".
- **The scoring components are not novel.** Accessibility–expression coactivity, a genomic
  distance prior, and TF motif evidence are, at the level of ingredients, what most of the
  methods below already combine. There is no formulation here that a reviewer would call new.

Bibliographic references are deliberately not invented; the one citation below was read
before citing.

---

## Summary

| Method | Unit of inference | Overlap with this work | Distance from this work |
|---|---|---|---|
| Signac `LinkPeaks` | peak–gene link | **Direct dependency** — supplies candidates and baseline | Zero. This repository is downstream of it |
| ArchR `Peak2GeneLinks` | peak–gene link | Same task, different implementation | Very close. Attempted here, excluded |
| Cicero | peak–peak co-accessibility | Adjacent; gene links are derived | Close |
| SCENT | peak–gene link, with inference | **Direct dependency** — external comparator | Zero. Used as the validator |
| SCARlink | gene expression from peaks | Same data, different objective | Close |
| CREMA | TF–site–gene circuit | Shares the TF/motif+coactivity ingredients | Close, different unit |
| TRIPOD | peak–TF–gene trio | Shares the trio structure | Moderate |
| SCENIC+ | eRegulon (TF + region + gene) | Shares ingredients, larger scope | Moderate |
| Pando | TF–peak–gene regulatory network | Shares ingredients | Moderate |
| LINGER | gene regulatory network, learned | Same data, much larger model | Further |
| FigR | DORC + TF regulator | Shares coactivity + TF logic | Close |

---

## Signac `LinkPeaks`

**What it does.** Correlates peak accessibility with gene expression across cells within a
distance window, calibrated against a background of peaks matched on GC content, accessibility
and width, returning a score, z-score and p-value per peak–gene pair.

**Why it matters.** It is the default in the Seurat/Signac ecosystem and therefore the most
widely used peak–gene linker in practice. It is also this repository's candidate generator and
baseline.

**How close.** Not a competitor. A dependency. Every candidate evaluated here came from
`LinkPeaks()` at a 500 kb window, and the `linkpeaks` score mode *is* its ordering.

**Claim to avoid.** Never "outperforms LinkPeaks" without qualification. The measured result is
that a reranking of LinkPeaks' own top 5,000 candidates concentrates SCENT-supported links higher
than LinkPeaks' ordering at top-100 and top-200. Recall is bounded by LinkPeaks; no link outside
its output was ever considered. Also avoid implying the comparison is symmetric — the baseline
chose the candidates, which cuts both ways and should be stated rather than glossed.

**Defensible positioning.** "A controlled reranking study over a LinkPeaks candidate universe,
with a distance-only control and distance-matched validation." The controlled part is the
contribution.

---

## ArchR `Peak2GeneLinks`

**What it does.** Correlates peak accessibility with gene expression across low-overlapping
aggregates of cells within a cis window, on ArchR's own Arrow-file infrastructure.

**Why it matters.** The main alternative to LinkPeaks, and the standard in ArchR-based workflows.
Its aggregate-based approach directly addresses the sparsity problem that limits per-cell
coactivity here.

**How close.** Very close — same task, same data type. ArchR was attempted as a benchmark
comparator during this project but did not produce benchmarkable output in this environment. It
is absent from the results.

**Claim to avoid.** Do not claim any comparison against ArchR. Do not imply its absence reflects
on ArchR — it reflects an integration failure in this environment. And do not claim the
aggregate-based approach was evaluated and rejected; it was not evaluated.

**Defensible positioning.** State plainly that ArchR was attempted, did not integrate stably, and
was excluded. If the standalone phase proceeds, ArchR aggregation is a directly relevant prior art
for the metacell work in `docs/future_standalone_v0.md` §4.

---

## Cicero

**What it does.** Estimates co-accessibility between pairs of ATAC peaks using aggregated
pseudo-cells and a graphical-lasso-style model, producing cis-co-accessibility networks. Gene
links are derived by connecting promoter peaks to co-accessible distal peaks.

**Why it matters.** One of the earliest and most influential single-cell cis-regulatory methods,
and the origin of the peak-aggregation idea now used widely.

**How close.** Close in spirit, different primary unit — peak–peak rather than peak–gene, and it
does not require paired RNA.

**Claim to avoid.** Do not claim to supersede Cicero. Do not claim the peak–peak structure was
tested; it was not. No Cicero comparison exists in this repository.

**Defensible positioning.** Different unit of inference. If Cicero is discussed at all, it should
be as prior art for aggregation, not as a benchmarked competitor.

---

## SCENT

**What it does.** Fits a Poisson regression of gene counts on peak accessibility per peak–gene
pair, with bootstrap inference, yielding an effect size and calibrated p-value. Designed for
single-cell multiome and intended to be robust to sparsity.

**Why it matters.** It provides *statistical* inference at the level of individual links, which
correlation-based methods do not. That is why it was chosen as the comparator here.

**How close.** Not a competitor. A dependency, used as the external comparator throughout
`docs/results_report.md`  — 52,482 tested rows, ~4,750 supporting.

**Claim to avoid.** Several, and these are the most important in this document.

- Do not claim to replace SCENT. This repository has no inferential machinery at all; it produces
  ordinal scores with no p-values or effect sizes.
- **Do not treat SCENT agreement as ground truth.** SCENT consumes the same RNA and ATAC matrices
  as the reranker. Agreement between two correlational methods on shared input is substantially
  weaker evidence than agreement with an orthogonal assay.
- Do not report SCENT support beyond its tested window. The sweep used 100 kb while candidates
  extend to 500 kb, so the `200_500kb` and `gt500kb` distance bins contain zero supported links
  for every method and their reported odds ratios are artifacts.
- Do not claim cell-type-specific concordance. SCENT ran with `scoring_celltype: ""`, i.e. a
  synthetic `all_cells` label.

**Defensible positioning.** "SCENT used as an independent correlational comparator and support
layer, with its window limitation stated." Nothing stronger.

---

## SCARlink

**What it does.** Predicts gene expression from tiled chromatin accessibility across a gene
locus using regularised regression, and derives per-region importance from the fitted model.

**Why it matters.** It reframes the problem as prediction rather than pairwise association,
which sidesteps the multiple-testing burden and captures joint contributions of multiple regions
to one gene.

**How close.** Close in data and objective, different in formulation. This repository scores
pairs independently; SCARlink models a locus jointly. It also produces a natural notion of
region importance, which is arguably a better-founded version of what a reranking score is
reaching for.

**Claim to avoid.** Do not claim competitiveness with regression-based prediction. No
predictive-accuracy evaluation exists here — there is no held-out expression prediction anywhere
in this repository. Do not present a multiplicative hand-set score as an alternative to a fitted
model.

**Defensible positioning.** Different formulation, and explicitly interpretable-by-construction
rather than fitted. That is a genuine trade-off, not an advantage, and should be described as
one.

---

## CREMA

**What it does.** Peak-agnostic. Scans the full TSS ±100 kb region for TF motif sites, measures
accessibility around each site, combines that with the TF's expression, and tests whether the
TF-expression × site-accessibility pattern associates with target gene expression. Output unit is
a TF–site–gene circuit.

**Why it matters.** It is the closest method in the same broad space — single-cell multiome
cis-regulatory inference using both modalities — and it uses the same three ingredients this
repository uses. Not being peak-dependent, it also avoids the peak-calling bottleneck.

**How close.** Close, but asking a different question:

```
CREMA:      Is TF X, acting through site Y near gene Z, associated with expression of Z?
this work:  Among candidate cis peak-gene pairs, which should be prioritized?
```

**Claim to avoid.** Do not present the TF/motif term as CREMA-like circuit inference. $T_p$
here is a **peak-level scalar with no gene dependence and no cell-type dependence**
(`docs/method_report.md` §8) — all pairs sharing a peak get the same value. It cannot express
TF-to-target specificity, which is precisely CREMA's object. Measured contribution is small and
inconsistent: `coactivity_tf` beats `coactivity` at top-100 (0.500 vs 0.460) but loses at
top-200 (0.515 vs 0.525) and in the 50–200 kb bin (1.930 vs 2.129).

**Defensible positioning.** The unit of inference differs — peak→gene prioritization versus
TF→site→gene circuit. If the standalone phase adds cell-type-aware TF evidence, keep the TF term
*supporting* a peak–gene link rather than becoming the object of inference, or the distinction
disappears.

---

## TRIPOD

**What it does.** Tests peak–TF–gene trios non-parametrically, using matched cell groups to
control for confounding, with explicit attention to marginal versus conditional associations.

**Why it matters.** It is unusually careful about confounding in trio inference, which is the
same class of problem this benchmark's proximity controls address.

**How close.** Moderate. Shares the trio structure but is statistically far more rigorous, and
its matching strategy is conceptually adjacent to this repository's distance-matched
stratification.

**Claim to avoid.** Do not claim comparable statistical rigour. The controls here are
descriptive — a distance-only ranking, binned odds ratios, threshold removal — not a formal
matched-inference framework.

**Defensible positioning.** Same confounding concern, addressed descriptively rather than
inferentially. TRIPOD's matching approach is relevant prior art for the standalone controls.

---

## SCENIC+

**What it does.** Joint analysis of scRNA and scATAC to build eRegulons — TF, enhancer region set,
and target gene set — via motif enrichment, region–gene links, and TF–gene correlation, with
downstream network and cell-state analysis.

**Why it matters.** The most complete pipeline in this space, widely adopted, and it delivers
interpretable regulon-level output that is directly usable biologically.

**How close.** Moderate. Region–gene linking is one internal component of SCENIC+, and its
ingredient list overlaps heavily with this repository's. Its scope is much larger.

**Claim to avoid.** Do not claim regulon or network inference. No TF-to-target-set inference
exists here, and no network is constructed at any point. Do not claim greater interpretability —
eRegulons are arguably more interpretable to a biologist than a multiplicative score.

**Defensible positioning.** This work addresses one narrow component that SCENIC+ contains, with
explicit proximity controls that pipeline-level tools do not typically report.

---

## Pando

**What it does.** Infers TF–peak–gene regulatory networks from multiome data using regularised
regression over TF-motif-containing peaks in a gene's cis region, on top of a Seurat-based
workflow.

**Why it matters.** Directly comparable ingredients — motif-containing peaks, cis window, joint
RNA/ATAC — with a fitted rather than hand-set combination, and it integrates cleanly into the
same ecosystem this repository uses.

**How close.** Moderate to close. Arguably the closest method to what the future standalone
direction would produce if it stayed on a ranking-accuracy axis.

**Claim to avoid.** Do not claim novelty for the combination of coactivity, motif evidence and a
cis window — Pando already fits that combination with learned weights. Do not claim an advantage
from hand-set parameters without evidence; $\lambda$ and $\alpha$ were not fitted or
held-out-selected.

**Defensible positioning.** Interpretable-by-construction with no fitted weights, benchmarked
against a distance-only control. The control, not the score, is the distinguishing feature.

---

## LINGER

**What it does.** Learns gene regulatory networks from single-cell multiome data using a model
that incorporates external bulk regulatory data, transferring learned regulatory relationships
to single-cell context.

**Why it matters.** It represents the learned, data-hungry end of the field and can leverage
information unavailable to a single dataset.

**How close.** Further away — much larger model class, external training data, network-level
output.

**Claim to avoid.** Do not position this repository against learned methods on accuracy. It has
five interpretable parameters, one dataset, and no training procedure. Do not claim that
interpretability compensates for that; it is a trade-off.

**Defensible positioning.** Deliberately at the opposite end of the complexity spectrum. A
transparent, fully inspectable score is a legitimate design point, and useful precisely because
its failure modes are visible — which is what let the proximity confound be detected here.

---

## FigR

**What it does.** Identifies domains of regulatory chromatin (DORCs) by finding genes with many
significantly correlated peaks, then associates TF regulators with DORCs using motif enrichment
and TF expression correlation.

**Why it matters.** Its DORC concept captures that regulatory signal is often distributed across
many peaks per gene rather than concentrated in one, and its TF-association step uses
essentially the same logic as this repository's TF term.

**How close.** Close on the TF/motif and coactivity logic; different in that it aggregates to
gene-level regulatory domains rather than scoring individual pairs.

**Claim to avoid.** Do not claim novelty for combining coactivity with motif-based TF
association — FigR does this. Do not claim per-link resolution is superior; DORC aggregation is a
deliberate noise-reduction choice, and given that sparsity is a known limitation here, it may be
the better choice.

**Defensible positioning.** Per-pair rather than per-domain resolution, with explicit distance
control. Note that FigR's aggregation logic is relevant prior art for the metacell work in
`docs/future_standalone_v0.md` §4.

---

## Benchmark precedent — BENGI

**What it is.** A Benchmark of candidate Enhancer-Gene Interactions, built by integrating the
ENCODE Registry of cCREs with experimentally derived interactions: RNAPII and CTCF ChIA-PET,
Hi-C, promoter-capture Hi-C, GEUVADIS and GTEx eQTLs, and crisprQTLs from a CRISPR screen —
21 datasets across 13 biosamples, over 162,000 cCRE–gene pairs, with negatives generated per
enhancer within a distance cutoff and cross-validation grouped by chromosome to prevent
overfitting.

**What it found.** A baseline distance method — rank pairs by inverse linear distance to the
gene's nearest TSS — outperformed both correlation-based unsupervised methods on every dataset
tested (average AUPR increase 127% over DNase-DNase, 77% over DNase-expression). The best
supervised method, TargetFinder, beat it within cell type but only modestly, and frequently
failed to beat it across cell types.

**Relation to this work.** Two distinct points.

*On the distance baseline.* The finding here that raw top-N support is proximity-driven is
consistent with BENGI's central result, which establishes it on experimental ground truth rather
than a correlational comparator. **This benchmark does not establish that distance is a strong
baseline; BENGI did.**

*On the setting.* BENGI's correlation methods correlate signals **across biosamples** — enhancer
accessibility against promoter accessibility, or accessibility against expression, over 32 and
112 cell types respectively. The paper attributes their failure to promoters being ubiquitously
active across cell types while enhancers are cell-type-specific, which washes out the
correlation. The coactivity score here is a **within-sample, across-cell** quantity, so that
mechanism does not transfer, and BENGI did not evaluate within-sample single-cell coactivity.
The within-bin distance-matched comparison in `docs/results_report.md` §5 is likewise not one
BENGI ran for its correlation methods — it applied distance-matched sampling to the supervised
model, where AUPR fell from 0.86 to 0.74 distance-matched and 0.61 promoter-distance-matched.

**Claim to avoid.** Do not claim the distance-baseline observation as novel. Do not claim
distance-matched stratification as a novel control — BENGI's Methods define five distance
quantiles and sample equally from each, the same move on a different substrate. And do not claim
this work contradicts BENGI; the correlation axis differs.

**Relevance to the next phase.** BENGI is the kind of validator
`docs/future_standalone_v0.md` §8 argues is required — experimental rather than correlational,
and not derived from the same two matrices as the score. Its crisprQTL set (K562, 4,937 tested
enhancers overlapping cCREs-ELS) is a distance-decoupled truth set that already exists in
curated form. Datasets and scripts: https://github.com/weng-lab/BENGI

Moore JE, Pratt HE, Purcaro MJ, Weng Z. *A curated benchmark of enhancer-gene interactions for
evaluating enhancer-target gene prediction methods.* Genome Biology 2020;21:17.
doi:10.1186/s13059-019-1924-8

## What positioning survives

Three statements are defensible on the current evidence.

1. **A controlled reranking benchmark with an explicit proximity control.** Ten score modes over
   one byte-identical 5,000-pair candidate universe, with a distance-only ranking, within-bin
   distance-matched enrichment, and proximal-removal thresholds. The control earns its place: it
   wins at top-50 on the unrestricted universe. Proximity controls are not novel — see the BENGI
   section above — but they remain uncommon in single-cell multiome reranking comparisons.

2. **A documented negative result about evaluation methodology.** Apparent improvements in
   peak–gene ranking can be produced by proximity collapse; SCENT-based validation is
   promoter-biased and window-limited; and a naive implementation of distance-matched enrichment
   emits meaningless odds ratios on empty bins. All three are demonstrated here with committed
   numbers. The distance-baseline component of this was established earlier and on stronger
   ground truth by BENGI; what is added here is the single-cell multiome setting and the
   window-alignment point.

3. **Reproducible infrastructure.** A containerised, `renv`-pinned, config-driven ablation
   harness with a chromosome-sharded SCENT sweep that completed genome-wide where a monolithic
   run was intractable. The harness is reusable independently of the scientific result.

## What positioning does not survive

- Any claim to be a peak–gene linking method.
- Any claim of novelty for the score family or its components.
- Any claim about distal enhancer–gene links — unmeasured, since the validator window is 100 kb.
- Any claim of superiority over the eleven methods above; only `linkpeaks` was benchmarked, and
  it also defined the candidate universe.
- Any claim of cell-type specificity.
- Any claim of validation against ground truth.
- Any claim that the proximity confound was eliminated. It was detected and quantified, which is
  the useful part, but `distance_only` remains the strongest method at the 50 kb threshold.

## Note on the future direction

If the standalone phase proceeds, the axis matters more than the implementation. Competing on
ranking accuracy puts it against eleven established methods with a simpler feature set and no
identified advantage. The alternative axis — cell-type-specific, tiered, **calibrated** output for
experimental follow-up — changes the required claim from "ranks better", which this evidence
cannot support, to "tier 1 has measurably higher orthogonal support than tier 3", which is
testable with the existing harness and does not require beating anyone. Most tools in the table
above output score tables; few output a tier with a stated hit rate and an explicit proximity
warning. See `docs/future_standalone_v0.md` §7.
