
## PBMC LinkPeaks reranker benchmark — current clean run

### Run status

Feature generation completed successfully after pinning local JASPAR2022 SQLite.

Candidate generation:
- LinkPeaks positive candidates after filtering: 15,806
- Feature-table candidates used for reranking: 5,000
- Candidate genes: 1,390
- Evaluation universe: same 5,000 peak-gene pairs for every score mode

Score modes run:
- linkpeaks
- coactivity
- distance_only
- coactivity_distance
- coactivity_tf
- full
- full_lambda_0_1
- full_lambda_0_2

### Main interpretation

The scoring idea has real signal.

Coactivity alone is strongly aligned with LinkPeaks and is not explained by distance-only ranking. Distance-only behaves like a promoter-proximity control and has very low top-link overlap with LinkPeaks.

TF/motif support changes rankings but does not dominate.

Distance prior is strong. `lambda_distance = 0.3` is probably too aggressive. `lambda_distance = 0.1` is the most reasonable current full-score candidate.

### Practical conclusion

Continue toward standalone cis peak-gene linker.

Recommended primary candidates:
- `coactivity_tf`: coactivity × TF support, no distance prior
- `full_lambda_0_1`: coactivity × mild distance prior × TF support

Sensitivity/control:
- `full_lambda_0_2`
- `full` / lambda 0.3
- `distance_only` negative/control

Do not claim `full` lambda 0.3 is best. Treat it as an aggressive distance-prior setting.

### Next step

Run SCENT as an external comparator on PBMC after committing this benchmark state.

Then run BC as second-dataset validation.

## Reranker biological validation plan

### Current status

The PBMC LinkPeaks reranker benchmark now produces a clean, reproducible result on a fixed candidate universe.

Current benchmark setup:

* LinkPeaks positive candidates after filtering: 15,806
* Feature-table candidates used for reranking: 5,000
* Candidate genes: 1,390
* Evaluation universe: same 5,000 peak-gene pairs for every score mode

Score modes already run:

* `linkpeaks`
* `coactivity`
* `distance_only`
* `coactivity_distance`
* `coactivity_tf`
* `full`
* `full_lambda_0_1`
* `full_lambda_0_2`

The current reranker should not be treated as the final method. It should be treated as a diagnostic scaffold for deciding whether the scoring components are worth carrying into a standalone cis peak-gene linker.

### Current interpretation

The reranker results support the following conclusions:

* Coactivity carries real signal.
* Distance-only is not sufficient and behaves mainly as a promoter-proximity control.
* TF/motif support changes rankings but does not dominate.
* Distance prior is strong and can easily overcompress rankings toward promoter-proximal links.
* `lambda_distance = 0.3` is likely too aggressive.
* `lambda_distance = 0.1` is the most defensible current full-score setting.
* Gene-only ORA is too weak to serve as the primary biological validation.

Therefore, the current result should be interpreted as:

> The scoring idea is biologically plausible and worth testing further, but the reranker has not yet proven superiority over LinkPeaks. The next validation must operate at the peak-gene link level, not only at the top-gene ORA level.

### Why gene-only ORA is insufficient

Gene-only ORA asks whether top-ranked genes are biologically plausible.

That is too weak because:

* PBMC top genes often produce generic immune terms regardless of link quality.
* ORA ignores which peak was linked to each gene.
* ORA discards peak-gene pair specificity.
* ORA is sensitive to the selected gene background.
* Distance/proximal bias can make gene-level enrichment look better without improving regulatory link quality.

Keep ORA only as a descriptive sanity check.

Do not use ORA alone to claim that the reranker is better.

### Main biological validation question

The correct validation question is:

> Are top-ranked peak-gene pairs biologically more plausible than lower-ranked candidate pairs and stronger than simple baselines?

The validation should compare:

* LinkPeaks ranking
* coactivity-only ranking
* coactivity + TF ranking
* mild full score: `full_lambda_0_1`
* stronger distance settings: `full_lambda_0_2`, `full`
* distance-only negative/control
* SCENT comparator output

### Validation 1: SCENT concordance

SCENT should be used as the primary external comparator.

Goal:

> Test whether reranked peak-gene links agree with an independent enhancer-gene method.

Do not use exact peak string matching as the main matching method, because SCENT and LinkPeaks/reranker peak sets may differ.

Use genomic overlap matching:

* represent peaks as genomic intervals
* match SCENT peaks to reranker peaks using `GRanges`
* use reciprocal overlap threshold, for example reciprocal overlap ≥ 0.5
* compare gene + overlapping-peak pairs

Metrics:

* SCENT-supported links in top 50 / 100 / 200
* fraction of top-ranked links supported by SCENT
* median rank of SCENT-supported links
* enrichment of SCENT-supported links in top decile vs all candidates
* AUPRC / AUROC if a clean binary SCENT-supported label can be constructed

Important: SCENT should not be run again as one monolithic all-genome job. Previous notes indicate that the all-cells native candidate run was too large. Use a tractable SCENT strategy:

* per-chromosome runs
* smaller cis windows
* stronger candidate prefiltering
* smoke-test scale before full scale

### Validation 2: distance-matched SCENT enrichment

Distance is a major confounder.

The key anti-bullshit test is:

> Within the same distance bin, are high-ranked reranker links more SCENT-supported than low-ranked links?

Distance bins:

* 0–10 kb
* 10–50 kb
* 50–200 kb
* 200–500 kb

For each bin, compare:

* top-ranked links by each score mode
* lower-ranked candidate links
* distance-only ranking
* LinkPeaks ranking

Useful metric:

* odds ratio for SCENT support in high-ranked vs background links within the same distance bin

This prevents the false conclusion that the model is biologically better when it is only choosing promoter-proximal links.

### Validation 3: promoted and demoted link inspection

Compare `full_lambda_0_1` against LinkPeaks.

Inspect:

* top links promoted by the reranker
* top links demoted by the reranker

For each inspected link, report:

* peak
* gene
* genomic distance
* LinkPeaks score
* coactivity score
* distance score
* TF/motif score
* top motif names
* whether the gene has plausible PBMC/immune relevance
* whether the peak has plausible regulatory support

This is not the main statistic, but it helps detect whether the reranker is doing biologically sensible prioritization or just moving arbitrary links.

### Validation 4: TF/motif support enrichment

The TF term should be validated directly.

Question:

> Do top reranked links have stronger motif/TF support than matched background links?

Compare top-ranked vs background candidate links for:

* `tf_score`
* `peak_motif_score`
* `peak_tf_score`
* top motif names
* motif support in distal links specifically

Important comparisons:

* `coactivity` vs `coactivity_tf`
* `coactivity_tf` vs `full_lambda_0_1`
* high-ranked distal links vs low-ranked distal links

If TF support does not improve link-level evidence or produces unstable rankings, keep it optional rather than central.

### Validation 5: distance and diversity behavior

For every score mode, report:

* top50 median distance
* top100 median distance
* fraction of distal links >50 kb
* unique genes in top100
* maximum number of top links assigned to one gene
* top-N overlap with LinkPeaks
* top-N overlap with SCENT-supported links if available

Interpretation:

* distance-only should behave as a negative/control ranking
* strong full models should not collapse completely to promoter-proximal links
* mild distance should preserve some distal regulatory candidates
* reranker should not simply concentrate all top links around a few genes

### Stop rule for reranker work

Do not keep tuning the reranker indefinitely.

Stop reranker polishing after the following validation package exists:

1. SCENT concordance
2. distance-matched SCENT enrichment
3. promoted/demoted link inspection
4. TF/motif support comparison
5. distance/distal/gene-diversity summary

Then make a decision.

### Decision rule

Continue toward standalone cis-window model if:

* `coactivity_tf` or `full_lambda_0_1` ranks SCENT-supported links above unsupported links
* it beats distance-only clearly
* it is not only winning in the promoter-proximal bin
* promoted links look biologically plausible
* TF/motif support improves or at least stabilizes biologically meaningful links

Do not continue polishing the reranker if:

* SCENT-supported links are not enriched in reranked top links
* distance-only explains most of the apparent signal
* improvements disappear after distance matching
* TF/motif support only adds noise
* promoted links look biologically arbitrary

### Project-level conclusion

The reranker is not the final method.

Its purpose is to identify which scoring components should be carried forward.

Current provisional design choices:

* keep coactivity as the core signal
* keep TF/motif support as a modifier
* keep distance as a mild prior only
* avoid aggressive distance weighting
* do not use ORA as the primary validation criterion

The next real method step, after the link-level validation package, is to replace LinkPeaks candidate generation with an own cis-window candidate generator and test whether the same scoring logic still works without LinkPeaks defining the candidate universe.

## Negative interpretation of current reranker benchmark

### Current conclusion

The current LinkPeaks-candidate reranker does **not** provide convincing evidence of improvement over LinkPeaks.

The reranker produces interpretable score behavior, but the biological validation does not support a strong positive claim. The most accurate interpretation is:

> The model reshuffles LinkPeaks-derived candidate links using biologically motivated features, but the reshuffling is not clearly beneficial.

This means the reranker should not be treated as the final method.

### What failed to look convincing

The main issue is that the modes expected to be most interesting do not clearly outperform the LinkPeaks baseline.

Observed pattern:

* `linkpeaks` remains biologically coherent by the available weak metrics.
* `coactivity` stays close to LinkPeaks but does not clearly improve it.
* `coactivity_tf` adds TF/motif influence but appears to weaken ORA-style biological signal.
* `full_lambda_0_1` gives a milder distance prior but still does not show clear improvement.
* `full_lambda_0_2` and `full` recover some gene-level enrichment but mainly by becoming more promoter-proximal.
* `distance_only` behaves as expected: mostly promoter-proximal and biologically weak.

The result is not random chaos, but it is also not a clear win.

### Main failure mode

The reranker appears to be doing biologically interpretable reshuffling rather than demonstrably better prioritization.

The score components behave in the intended direction:

* coactivity contributes real ranking structure
* TF/motif support shifts rankings
* distance regularization strongly affects prioritization
* distance-only is a poor control

However, the combined score does not produce a convincing biological improvement over LinkPeaks.

This suggests that the current reranking formulation is not strong enough as a standalone result.

### Why this is not publishable as-is

The current result is too easy to criticize.

Main weaknesses:

* No clear improvement over LinkPeaks.
* Gene-only ORA is weak and cannot validate peak-gene links.
* Top-ranked reranker outputs may lose immune/PBMC enrichment compared with LinkPeaks.
* Stronger distance settings improve some summaries but mostly by pushing links closer to promoters.
* LinkPeaks defines the candidate universe, so the reranker is constrained by and compared back to the same method.
* The most promising modes are not obviously better than the baseline.
* The result supports “interpretable behavior,” not “better method.”

Therefore, the reranker should not be framed as a successful method.

### What the benchmark still taught

The work was not completely useless.

It established several important points:

* The pipeline now runs reproducibly.
* The candidate universe is controlled.
* Score modes can be compared on the same 5,000 peak-gene pairs.
* Distance-only is a useful negative/control ranking.
* Coactivity is a real signal, but not enough to prove improvement.
* TF/motif support can modify rankings, but does not rescue the current model.
* Distance must be used carefully; aggressive distance weighting collapses toward promoter-proximal links.
* ORA is not a sufficient validation strategy for peak-gene linking.

This is a useful diagnostic result, but not a successful reranker result.

### Decision

Stop polishing this reranker.

Do not spend more time on:

* more lambda tuning
* more alpha tuning
* more ORA interpretation
* trying to force the reranker to look better
* treating this as the final method

The reranker should be frozen as a prototype / diagnostic branch.

### Honest project status

Current status:

> The LinkPeaks reranker did not demonstrate a convincing biological improvement. It should be treated as a negative or weak diagnostic result. The useful pieces are the feature calculations, evaluation harness, and lessons about score behavior, not the reranker itself.

### Possible next paths

There are only three rational next choices.

#### Option A — Archive / stop

Stop the project here.

This is defensible if the goal was specifically to produce a LinkPeaks reranker.

Conclusion:

> Reranking LinkPeaks candidates with the current coactivity × distance × TF/motif score was not strong enough to justify further development.

#### Option B — Audit for bugs first

Before fully trusting the negative result, run a systematic audit.

Check:

* candidate pairs are identical across modes
* score formulas are implemented exactly
* rank columns match score order
* RNA and ATAC cell orders are aligned
* peak identifiers match TF/motif scores
* distance scores decrease correctly with genomic distance
* no stale outputs are being interpreted

If a systematic bug is found, fix and rerun.

If no bug is found, accept the weak result.

#### Option C — Stop reranker work and test standalone candidate generation

Only continue if the goal changes from “rerank LinkPeaks” to “build an own cis-window peak-gene linker.”

The next valid method question would be:

> Does the same coactivity / TF / mild-distance logic work when LinkPeaks no longer defines the candidate universe?

This would require a separate minimal standalone candidate-mode test.

But this should not be framed as rescuing the reranker. It is a different method direction.

### Final note

The current reranker result is not encouraging.

The strongest honest conclusion is:

> The current reranker produces interpretable but unconvincing reshuffling of LinkPeaks candidates. It should not be optimized further as the main method. Future work should either audit for systematic bugs, stop the project, or move to a genuinely standalone candidate-generation model.

## Modified distance ablation finding

### Purpose

A modified distance prior was added to test whether the distance component could be made less aggressively promoter-biased while still contributing useful genomic plausibility.

The original distance formulation was:

```text
old_distance_modifier = (1 - lambda_distance) + lambda_distance * distance_score
```

This only penalizes distal links. Nearby links are preserved, while distant links are downweighted.

The modified distance formulation was:

```text
modified_distance_modifier = 1 + lambda_distance * (distance_score - 0.5)
```

This makes distance act more symmetrically:

```text
nearby links  -> mild boost
mid-distance  -> approximately unchanged
distal links  -> mild penalty
```

The goal was not to make distance dominate, but to test whether a smoother distance prior could retain more distal regulatory candidates than the original full model.

### New ablations added

The following new modes were added:

```text
distance_mod_only_lambda_0_1
full_moddist_lambda_0_1
full_moddist_lambda_0_2
```

The main intended comparison set is now:

```text
linkpeaks
coactivity
coactivity_tf
full_lambda_0_1
full_lambda_0_2
full_moddist_lambda_0_1
full_moddist_lambda_0_2
distance_only
distance_mod_only_lambda_0_1
```

### Main finding

The modified distance prior behaves correctly, but it does not radically change the ranking at lambda 0.1.

`full_moddist_lambda_0_1` is very similar to `full_lambda_0_1` in top-ranked behavior. This is expected because the near-vs-far ranking pressure is almost the same at lambda 0.1 under the old and modified formulas.

At lambda 0.2, the modified distance prior is somewhat less harsh than the original distance prior. It preserves slightly more distal links than `full_lambda_0_2`, but it still shifts rankings strongly toward promoter-proximal links.

### Distance behavior summary

Observed pattern:

```text
coactivity_tf:
  no distance prior
  preserves the most distal links among plausible reranker variants

full_lambda_0_1:
  mild original distance prior
  reduces distal fraction noticeably

full_moddist_lambda_0_1:
  modified distance prior
  very similar to full_lambda_0_1

full_lambda_0_2:
  stronger original distance prior
  strongly promoter-proximal

full_moddist_lambda_0_2:
  less harsh than full_lambda_0_2
  but still substantially distance-biased

full / lambda 0.3:
  too aggressive
  should be treated as a distance-pressure sensitivity run, not as the preferred model
```

Approximate top-ranked behavior from the PBMC validation run:

```text
coactivity_tf:
  top50 median distance ~11.3 kb
  top50 distal fraction >50 kb ~32%

full_lambda_0_1:
  top50 median distance ~5.9 kb
  top50 distal fraction >50 kb ~20%

full_moddist_lambda_0_1:
  top50 median distance ~5.9 kb
  top50 distal fraction >50 kb ~20%

full_lambda_0_2:
  top50 median distance ~3.1 kb
  top50 distal fraction >50 kb ~6%

full_moddist_lambda_0_2:
  top50 median distance ~3.8 kb
  top50 distal fraction >50 kb ~10%
```

### Interpretation

The modified distance prior does not invalidate the earlier conclusion that distance is a strong and potentially dangerous component.

Distance is clearly doing real work. It changes rankings and pushes top links toward shorter genomic distances. This is biologically plausible up to a point, but it can easily become promoter-proximity bias.

The motif/TF component is also doing real work. `coactivity_tf` changes the ranking relative to coactivity alone and promotes links with higher TF/motif support. This supports keeping motif/TF support as a real scoring component rather than treating it as cosmetic.

The current evidence suggests:

```text
coactivity = core signal
TF/motif support = useful modifier
distance = useful but risky prior
modified distance = cleaner formulation, but not a decisive improvement yet
```

### Current method decision

For the next validation step, the primary full-method candidate should be:

```text
full_moddist_lambda_0_1
```

Reason:

```text
It matches the intended biological model:
coactivity × TF/motif support × mild distance prior
```

The most important ablation should be:

```text
coactivity_tf
```

Reason:

```text
It tests whether adding distance improves or hurts after coactivity and TF/motif support are already present.
```

The required controls are:

```text
linkpeaks
coactivity
distance_only
```

Reason:

```text
linkpeaks = source-method baseline
coactivity = core signal only
distance_only = promoter-proximity control
```

### What not to claim

Do not claim that the modified distance model has proven superiority over LinkPeaks.

Do not claim that the distance prior is solved.

Do not claim that lambda 0.2 or lambda 0.3 are preferred.

Do not use ORA alone as evidence that distance improves the model.

### SCENT decision rule

The next validation should ask:

```text
Does full_moddist_lambda_0_1 rank SCENT-supported links above unsupported links?
```

Critical comparison:

```text
full_moddist_lambda_0_1 vs coactivity_tf
```

Interpretation:

```text
If full_moddist_lambda_0_1 beats coactivity_tf after distance matching:
  the mild distance prior is adding useful biological signal.

If coactivity_tf beats full_moddist_lambda_0_1 after distance matching:
  distance is mostly adding bias.

If both beat distance_only:
  coactivity/TF signal is likely real.

If neither beats LinkPeaks or distance_only:
  stop tuning the reranker and salvage the feature engineering only.
```

### Practical conclusion

The modified distance experiment was useful. It showed that changing the distance formulation can reduce harshness at stronger lambda values, but lambda remains the main sensitivity point.

The best current full-method candidate is:

```text
full_moddist_lambda_0_1
```

The best no-distance comparator is:

```text
coactivity_tf
```

The next step is not more distance tuning. The next step is SCENT validation with distance-matched analysis.

## PBMC SCENT validation summary

SCENT validation was run genome-wide across chr1–chr22 using the completed SCENT chromosome sweep. The sweep produced 52,482 tested SCENT peak–gene rows, of which 4,758 passed the positive/significant support rule (`pvalue_positive`). Ranked methods were restricted to the chromosomes covered by the SCENT sweep before comparison.

The main result is that the full reranker models outperform the original LinkPeaks ranking by SCENT link-level support. At top 200 links, SCENT-supported fractions were:

* `full_lambda_0_1`: 0.605
* `full_moddist_lambda_0_1`: 0.600
* `coactivity`: 0.525
* `coactivity_tf`: 0.515
* `distance_only`: 0.510
* `linkpeaks`: 0.445

At top 100, the same pattern mostly holds:

* `full_lambda_0_1`: 0.580
* `full_moddist_lambda_0_1`: 0.570
* `distance_only`: 0.530
* `coactivity_tf`: 0.500
* `coactivity`: 0.460
* `linkpeaks`: 0.430

At top 50, `distance_only` was highest, but this is not convincing biological evidence because the top distance-only links are almost entirely promoter/TSS-proximal. Its median distances were extremely small: about 3.5 bp at top 50, 7.5 bp at top 100, and 15 bp at top 200, with promoter fraction near 1.0. This shows that SCENT support is strongly enriched near promoters, but also that a pure distance model collapses onto trivial proximal links.

The most important result is the distance-matched enrichment analysis. Within distance bins, the full models still enrich for SCENT-supported links over lower-ranked links. In the 0–10 kb bin, LinkPeaks had an odds ratio of about 2.60, while coactivity, coactivity_tf, full_lambda_0_1, and full_moddist_lambda_0_1 were around 4.68–5.06. In the 10–50 kb bin, LinkPeaks had an odds ratio of about 1.77, while the full models were about 3.15. In the 50–200 kb bin, LinkPeaks was about 1.65 and the full models were about 2.23. This suggests the reranker is not only exploiting global proximity; its coactivity/TF/full scores enrich SCENT-supported links even after distance stratification.

The comparison between `full_lambda_0_1` and `full_moddist_lambda_0_1` does not show a meaningful empirical difference. Their top-200 pair overlap was 199/200, and their SCENT support metrics were nearly identical. Therefore, `full_moddist_lambda_0_1` should be kept only as the cleaner distance formulation, not because it clearly outperformed the original λ=0.1 distance model.

Overall, SCENT gives a more favorable result for the reranker than the earlier gene-level ORA analysis. ORA favored LinkPeaks and was weak/gene-only, whereas SCENT is link-level and supports the diagnostic value of the reranker components. The best interpretation is that coactivity, TF/motif support, and a mild distance prior all carry useful information, but distance alone is not a satisfactory model because it over-prioritizes promoter-proximal links.

Conclusion: the LinkPeaks-candidate reranker should still be treated as a diagnostic scaffold, not a final method. However, the SCENT validation supports carrying the main components forward into the standalone cis-window prototype: coactivity, TF/motif support, mild distance prior, and distance-only as a control. A promoter-proximal dampening or hump-shaped distance prior may be worth testing as an ablation in the standalone version.

### removing proximal links

After excluding links within 10 kb or 25 kb of the TSS, full_lambda_0_1 and full_moddist_lambda_0_1 still outperform LinkPeaks by SCENT-supported top-N fraction. This argues that the full-model advantage is not solely due to promoter/TSS-proximal links. However, after excluding links within 50 kb, the advantage over LinkPeaks becomes modest, and distance-only remains high by selecting links just above the cutoff. Thus, the current reranker signal appears strongest for proximal-to-intermediate cis links rather than clearly distal enhancer-gene links.