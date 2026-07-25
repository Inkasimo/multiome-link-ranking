# Assistant handoff summary

Decisions and context from the review session of 2026-07-25. **This is not a summary of the
repository** — for that read `README.md`, `docs/method_report.md` and `docs/results_report.md`.
The purpose here is to preserve settled decisions so they are not relitigated.

If you are an assistant picking this up: the owner's goal is to **finish and freeze** this repo,
not to improve the science. Resist scope expansion. Most open questions below are deliberately
deferred, not overlooked.

---

## 1. What the repo is

A PBMC LinkPeaks-candidate peak–gene **reranking benchmark**. Ten interpretable score modes
ranked over one fixed 5,000-pair candidate universe (4,976 after restriction to SCENT-covered
chromosomes), with SCENT as an external comparator and explicit proximity controls.

Dataset: 10x `pbmc_unsorted_10k`, Cell Ranger ARC 2.0.0, `GRCh38-2020-A`, CC BY 4.0.

## 2. What it is not

Not a peak–gene linking method. Not causal inference. Not an enhancer–gene atlas. Not
cell-type-specific. Not a replacement for Signac, ArchR, Cicero, SCENT, SCARlink, CREMA, TRIPOD,
SCENIC+, Pando, LINGER or FigR. Recall is bounded by LinkPeaks, which both generated the
candidates and supplies the baseline.

## 3. Safest framing — settled

"Benchmark", never "method". The word appears in the repo name, the `CITATION.cff` abstract, and
the README's first line. Never "outperforms" unqualified, never "validated".

## 4. Main value proposition — settled

In order: (a) the proximity control suite — `distance_only` as an explicit mode, distance-matched
stratification, proximal removal; (b) the documented negative result about proximity confounding
in peak–gene evaluation; (c) the reproducible harness, including chromosome-sharded SCENT.

**The score itself is not the value proposition.** Its components are the field's consensus
ingredients and it does not cleanly beat a distance prior.

## 5. Proximity interpretation — settled position

- Distance-matched enrichment is the real finding: 5.06 / 3.15 / 2.23 against LinkPeaks'
  2.60 / 1.78 / 1.65, with `distance_only` at 1.59 / 0.92.
- The raw top-N advantage is **partly** a proximity effect — the full models' top-100 median
  distance is 7.3 kb against LinkPeaks' 16.2 kb.
- `distance_only` wins outright at top-50 (0.580 vs 0.440, median distance 3.5 bp) and beats every
  model at the 50 kb removal threshold. **The confound is displaced by the controls, not
  eliminated.**
- The 0–10 kb enrichment is a **coactivity + TF** result. `coactivity_tf`, `full_lambda_0_1` and
  `full_moddist_lambda_0_1` are all 5.057 there because \(f_\lambda\) is constant when
  \(D \approx 1\). The distance prior contributes nothing to the headline.
- SCENT tested only ±100 kb while candidates reach 500 kb, so the distal regime is **unmeasured**.
  The `200_500kb` and `gt500kb` odds ratios (8.906, 0.333) are empty-cell artifacts.
- Gene-level ORA favours LinkPeaks, 17 terms to 5. Superseded by the link-level SCENT analysis,
  **not refuted**. Both stay in the record.

## 6. Candidate universes — the distinction that keeps getting lost

| Set | Size | Note |
|---|---|---|
| LinkPeaks candidates | 15,806 | after `--candidate-filter positive_score` |
| Feature table / reranked | 5,000 | top-k by `link_score` |
| SCENT-evaluated | 4,976 | after `restrict_to_scent_chrs: true` |
| De novo cis-window | 117,811 | `run_scent_chr_sweep.R:325–365`, own code, 36.5% overlap with LinkPeaks |

Never mix 5,000 and 4,976. **Peak ID formats differ on disk** — `chr1-x-y` in LinkPeaks-derived
files, `chr1:x-y` in SCENT candidate files. Joining on the raw string returns zero matches. Split
on `[:-]` first.

## 7. Doc corrections outstanding

See `docs/DOC_CORRECTIONS.md` — five edits, all measured, none needing a re-run. Apply and stop.
Also: `README.md` links to `docs/lab_notebook.md`, which only exists after the
`developer_notes.md` rename.

## 8. Cleanup order

`TODO.md` §5. **`TODO.md` §0.1 gates everything** — the Snakefile calls
`scripts/run_linkpeaks_reranker.R`, which does not exist. The rename from
`run_linkpeaks_reranker_without_scent.R` is **decided and verified** by CLI-contract diff: all 19
Snakefile flags accepted, output filenames match lines 133–136, and neither script contains any
SCENT code (the suffix is a misnomer). Nothing can run until this is done.

## 9. Claims that must be avoided

The 8.906 odds ratio (empty-cell artifact). Anything about distal links. "Validated" or "ground
truth" for SCENT — it is correlational, shares both input matrices, and is promoter-biased.
Cell-type specificity — SCENT ran with a synthetic `all_cells` label and no score is stratified.
Novelty for the score family. `full_moddist` as an empirical improvement — 199/200 top-200 overlap
with `full_lambda_0_1`, and the algebra predicts it (\(f^{mod} = f^{orig} + \lambda/2\)).
Absolute SCENT-supported fractions as hit rates — no multiplicity control across 52,482 tests.

## 10. What standalone v0 would test — future work only

**Not implemented. Do not start before the repo is frozen and released.**

The de novo cis-window candidate generator already exists (§6). The real objective is not a
twelfth linker competing on ranking accuracy — it is cell-type-specific **calibrated tiers** for
experimental follow-up, which changes the claim from "ranks better" (unsupportable) to "tier 1 has
higher orthogonal support than tier 3" (testable).

Two traps: the candidate set and the SCENT validator are both ±100 kb, so scoring one against the
other makes the confound worse with no proximal-removal escape route — break the symmetry
deliberately. And cell-type stratification needs annotation (currently unlabelled clusters) plus
metacell aggregation (coactivity gets noisier when cells are split), neither of which exists.

The deeper constraint, stated once: every validator available is correlational and derived from the
same two matrices. An orthogonal, distance-decoupled truth set is what would make any result
defensible. That likely means leaving PBMC.

## 11. Go / no-go

Full criteria in `docs/future_standalone_v0.md` §6. Summary: **go** only if tiers separate on a
self-generated candidate universe and the full model beats `distance_only` — not merely LinkPeaks —
within distance bins after 25 kb proximal removal. **No-go** if tiers do not separate globally;
they will not separate after splitting cells by type.

Time-box any continuation to two weeks with the stop rule written down first. A no-go outcome
still improves the public repo, because it closes the question the current benchmark leaves open.

## 12. Settled — do not relitigate

- Dataset is `pbmc_unsorted_10k`, recovered from the fragment-file header. **Not** the Signac
  vignette's `pbmc_granulocyte_sorted_10k`, despite identical local filenames. CC BY 4.0.
- Lead with SCENT; ORA is superseded but not refuted.
- Publish `docs/lab_notebook.md` with a header marking it unedited working notes.
- Script rename decided (§8).
- Benchmark-not-method framing decided (§3).
- **No new analysis before release.** Four checks were considered and deliberately deferred:
  activity-matched enrichment, `adj` score modes, re-evaluation on all 15,806, and `tf_score`
  versus peak width. Each is documented as untested. See `docs/DOC_CORRECTIONS.md` §"Not doing".
- `adj` is a **failed correction**, not an unused improvement. Spearman −0.867 with marginals and
  −0.735 with its own numerator `mul_strict`; ratio normalisation against a 0.968-collinear
  denominator inverts the statistic and selects rare-feature pairs. Do not add it as a robustness
  ablation.

## 13. Known unknowns

JASPAR2022 sqlite redistribution terms (only remaining release blocker). Author ORCID and
affiliation. Final repo and image names — three names currently disagree
(`multiome-link-ranking` in the Dockerfile, `multiome-link-ranking:pilot` in `run_analysis.py`).
Cell count after QC — **there is no QC step**, so N = barcodes in the filtered matrix ≈ 12,012.
Total runtime — all sweep rows say `skipped_existing` with blank `runtime_minutes`. The
`--min-distances` / `--high-fraction` values behind the committed control outputs. Contents of
`results/pbmc/audits/reranker_output_audit_checks.csv`, excluded from the handoff tarball.

## 14. Two known defects, documented not fixed

**No cell-level QC.** `CreateSeuratObject()` at `run_linkpeaks_reranker.R:380` takes no
`min.cells`/`min.features` and there is no `subset()`. All ~12,012 barcodes are used. Fixing this
would invalidate every committed result — v0.2 only.

**Annotation mismatch.** Counts come from `GRCh38-2020-A`; TSS coordinates driving `distance_bp`
come from `EnsDb.Hsapiens.v86` (Ensembl 86, 2016). The Signac vignette pairs the same two sources,
so this is convention rather than error, but it is real coordinate drift on the central variable.
