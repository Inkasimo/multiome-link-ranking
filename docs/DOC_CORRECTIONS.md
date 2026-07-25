# Doc corrections — apply these five, then stop

Five edits. All measured, none requiring a re-run. Apply, commit, move on to cleanup.

If you have not yet edited the generated docs locally, you can skip this file and just take
the updated copies — these corrections are already applied there.

---

## 1. `docs/results_report.md` §5 — the tight-bin result is coactivity+TF, not the full model

**Why:** in the 0–10 kb bin `coactivity_tf`, `full_lambda_0_1` and `full_moddist_lambda_0_1` all
score **5.057**, identical to three decimals, because \(D \approx 1\) there so \(f_\lambda\) is
constant. The distance term contributes nothing to the strongest result in the study.

**Insert** after the "Interpretable bins" table, before the paragraph beginning "This is the
strongest evidence in the benchmark":

> Note that in the `0_10kb` bin, `coactivity_tf`, `full_lambda_0_1` and
> `full_moddist_lambda_0_1` are numerically identical at 5.057. Within that bin
> \(d \ll d_0\), so \(D \approx 1\) and the distance modifier is effectively constant across
> candidates — the three modes reduce to the same function. The tight-bin enrichment is
> therefore a **coactivity + TF/motif** result, not a result about the distance prior. The
> distance term only begins to differentiate in the wider bins (10–50 kb: 3.017 for
> `coactivity_tf` against 3.148 for the full models). Any statement that "the full model shows
> the strongest distance-matched enrichment" should be read with this in mind.

---

## 2. `docs/method_report.md` §6 — coactivity is associated with detection rate

**Insert** at the end of §6, after the paragraph beginning "Two properties of \(A_{pg}\)":

> **\(A_{pg}\) is associated with marginal detection rates.** Recovering the marginal activity
> product from the stored columns as `mul_strict / adj` gives
> Spearman(`mul_weigh`, marginals) = **+0.682** over the 5,000 candidates. `mul_weigh` has no
> natural null: under independence of two standardised variables its expectation is positive and
> depends on the detection rate of each feature, so pairs where both gene and peak are detected
> in an intermediate fraction of cells score higher for reasons unrelated to coactivity.
> LinkPeaks controls for exactly this through a background matched on accessibility and GC;
> `mul_weigh` does not. **This is not conditioned on anywhere in the benchmark.**
>
> The feature table contains an adjusted variant, `adj` = `mul_strict` divided by the marginal
> activity product — a lift ratio, null at 1. It is stored but **not used by any score mode**,
> and it was inspected in early prototype work without being adopted. Re-examination shows the
> ratio normalisation overcorrects rather than corrects: Spearman(`adj`, marginals) = **−0.867**
> and Spearman(`adj`, `mul_strict`) = **−0.735**, i.e. `adj` is negatively correlated with its
> own numerator. Dividing by a denominator that is 0.968-collinear with the numerator produces a
> ratio governed by the denominator, so ranking by `adj` preferentially selects rare gene × rare
> peak pairs. It is recorded here as a failed correction, not a robustness result.

**Add** to §14, in the "Methodological" block, renumbering the rest:

> **Coactivity is not conditioned on marginal activity.** `mul_weigh` correlates with the
> marginal detection-rate product at +0.682. Whether the coactivity term's contribution survives
> conditioning on marginal activity **is not tested in this release**. The natural check — an
> activity-matched enrichment analysis mirroring the distance-matched one in §12 — is deferred.

---

## 3. `docs/results_report.md` §11 — SCENT support has no multiplicity control

**Add** as a new numbered item in "What cannot be claimed":

> **Not** that the SCENT-supported labels are individually reliable. The support rule
> (`beta > 0` and `boot_p ≤ 0.05`) was applied to 52,482 tests with no correction for multiple
> testing. Against a one-sided 2.5% null expectation of roughly 1,300 rows, the observed 4,758
> supporting rows represent about 3.6-fold enrichment — but a substantial minority of individual
> support labels, plausibly a quarter to a third, may be false positives. Because all methods are
> scored against the same labels, the *comparisons* between methods remain valid and if anything
> are attenuated toward the null. The *absolute* fractions are not. "60.5% of the top 200 are
> SCENT-supported" must not be read as "60.5% are real links."

---

## 4. `docs/results_report.md` — add a new §8.1 for objections that were tested and came back clean

**Insert** as a subsection at the end of §8, "Distance-only confounding":

> ### 8.1 Two objections tested and not supported
>
> Both of these are obvious challenges to the distance-matched result. Both were checked against
> the committed feature table and neither holds. They are recorded because a negative check is
> worth as much as a positive one.
>
> **Coactivity does not smuggle in proximity.** If `mul_weigh` were itself distance-dependent,
> the distance-matched analysis would be circular. It is not:
> Spearman(`mul_weigh`, `distance_bp`) = **−0.079**, and the marginal detection-rate product is
> uncorrelated with distance at **−0.001**. Proximity enters the score only through the explicit
> distance term and through the composition of the candidate set — which is what binning on
> distance controls.
>
> **Truncating to the top 5,000 does not disadvantage the baseline.** `link_score` is heavily
> right-skewed, so removing the low-scoring bulk leaves dispersion essentially unchanged: SD
> 0.0425 in the top 5,000 against 0.0415 across all 15,806, a ratio of 1.025. A Thorndike
> range-restriction correction gives 0.98×, i.e. no attenuation. More directly, the universe was
> selected *by* `link_score` and the `linkpeaks` mode ranks *by* `link_score`, so its top-200 is
> identical to its top-200 over the full 15,806 — its ranking is invariant to the truncation and
> cannot be biased by it.

---

## 5. `README.md` — one line in Limitations

**Add** after the bullet beginning "No cell-type stratification anywhere":

> - The coactivity term `mul_weigh` is associated with marginal gene and peak detection rates
>   (Spearman +0.682) and is not conditioned on them. An adjusted variant `adj` exists in the
>   feature table but overcorrects and is unused. Whether the coactivity contribution survives
>   conditioning on marginal activity is not tested in this release.

---

## Not doing, deliberately

Recorded so it is not relitigated later. Each is legitimately deferrable; a benchmark is allowed
to state that something was not tested.

| Deferred | Why it is safe to defer |
|---|---|
| Activity-matched enrichment | Documented as untested in §14 and the README |
| `adj` / `adj_tf` / `adj_full_*` score modes | `adj` is a failed correction; adding it as a robustness ablation would misrepresent it |
| `mul_strict` as a score mode | Spearman +0.968 with marginals — it *is* a detection-rate statistic, nothing to learn |
| Re-evaluation on all 15,806 candidates | Range restriction shown negligible (SD ratio 1.025); this is a new experiment about complementarity, not a fix |
| `tf_score` vs peak-width check | Documented as an untested possibility in the method report |
| Residualised coactivity mode | v0.2 |

Anything on this list that later seems urgent: reread `docs/assistant_handoff_summary.md` §11
before starting.
