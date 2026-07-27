# Future standalone v0 — plan

**Nothing in this document is implemented in this repository.** The current repository is a
reranking benchmark over a LinkPeaks-defined candidate universe. The standalone method described
here is the broader scientific goal and the next phase of work. Keeping the two separate is the
point of documenting them separately.

This is deliberately shorter than `docs/method_report.md` and `docs/results_report.md`, because
those describe results and this describes intentions.

---

## 1. Why standalone is the next step

The current benchmark cannot answer its own most interesting question. LinkPeaks both defines
the candidate universe and supplies the baseline ranking, so:

- Recall is bounded by LinkPeaks at every score mode.
- "Beats LinkPeaks" means "reorders LinkPeaks' own output", not "finds links LinkPeaks missed".
- Whether the scoring carries signal when LinkPeaks is not choosing the candidates is
  structurally untestable here.

That last point is the scientific reason to go standalone. The practical goal, however, is not
another ranking method — see §7.

---

## 2. What the reranker taught

Six findings that constrain the next phase. Sources in `docs/results_report.md`.

1. **Coactivity is the dominant signal.** `coactivity` alone raises the top-200
   SCENT-supported fraction from 0.445 (LinkPeaks) to 0.525, and its within-bin odds ratios
   (4.68, 2.89, 3.08) are well above the baseline's (2.60, 1.78, 1.94). In the outermost
   testable bin, 50–100 kb, coactivity alone is the **strongest** method — ahead of every
   composite score.
2. **The distance prior is inert where it matters, and the wrong shape.** \(\lambda = 0.1\)
   and \(\lambda = 0.3\) give identical within-bin odds ratios in every proximal bin, and
   \(\lambda = 0.3\) is slightly worse at 25–50 kb and 50–100 kb. Raising \(\lambda\) buys
   raw support only by becoming more promoter-proximal. Separately, top-decile enrichment is
   **not monotone in distance** — 25–50 kb (4.21–4.66) exceeds 10–25 kb (2.20–2.55) — while the
   prior decreases monotonically. For v0 the defensible primary is \(\lambda = 0\), with
   distance retained purely as a stratifier and control, and \(\lambda = 0.1\) kept only for
   continuity with these results. If a distance term returns it should be the hump-shaped
   variant in §4, with \(d_0\) fitted against the fine bins rather than hand-set.
3. **Proximity dominates raw top-N selection but carries no within-bin signal.**
   `distance_only` wins at top-50 on the unrestricted 500 kb universe with a median distance of
   3.5 bp, so any support fraction quoted without a distance control is uninterpretable.
   Restricted to SCENT's 100 kb tested window, however, `distance_only` is the **weakest**
   method at every proximal-removal threshold and depth, and its distance-matched odds ratio is
   **below 1** in two fine bins (0.840 at 10–25 kb, 0.687 at 25–50 kb). Ranking by proximity
   inside a distance bin is worse than arbitrary. The control suite — distance-only ranking,
   distance-matched stratification, proximal removal, and restriction to the validator's tested
   window — is the most reusable output of this work and should be carried forward unchanged.
   *An earlier version of this section stated that `distance_only` beats every model after 50 kb
   removal. That was an artifact of scoring untested distal candidates as unsupported and is
   retracted; see `docs/results_report.md` §6.*
4. **Peak-level TF/motif support helps proximally and costs distally.** \(T_p\) has no gene and
   no cell-type dependence. `coactivity_tf` improves on `coactivity` at 0–10 kb (5.06 vs 4.68)
   and 10–25 kb (2.55 vs 2.20), but is worse at 50–100 kb (2.64 vs 3.08) and at top-200 on the
   raw universe. The term behaves as a promoter-context proxy rather than as evidence of
   TF-to-target regulation, and the sign flip is the clearest available argument that a gene-
   and cell-type-aware TF term is required for this component to be worth its complexity.
5. **The distance reparameterisation is empirically inert.** `full_moddist` overlaps
   `full_lambda_0_1` at 199/200 in the top 200, with identical odds ratios. Keep the cleaner
   form; do not spend more time on distance-prior shape without a better validator.
6. **The evaluation axis is the bottleneck, not the score.** Every result was limited by what
   SCENT could test, not by the scoring. More \(\lambda\) or \(\alpha\) tuning cannot improve
   the evidence.

---

## 3. Own cis-window candidate generation — largely already built

The candidate-generation milestone is closer to complete than the roadmap assumed.
`scripts/run_scent_chr_sweep.R` lines 325–365 already implement the full specification, in this
repository's own `data.table` code, with no LinkPeaks or SCENT involvement:

| Requirement | Implementation |
|---|---|
| **Expressed genes** | `expr_frac_gene <- Matrix::rowMeans(rna_counts > 0)`, retain \(\geq\) `min_pair_frac` (`:325, :328`) |
| **Accessible peaks** | `expr_frac_peak <- Matrix::rowMeans(atac_counts > 0)`, retain \(\geq\) `min_pair_frac` (`:326, :329`) |
| **Same chromosome** | `data.table` merge on chromosome (`:353–358`) |
| **TSS ±100 kb window** | `distance_bp := abs(peak_mid - tss)`, retain \(\leq\) `link_distance` (`:364–365`) |

Its output is committed as 22 files, `results/pbmc/scent_chr_sweep_100kb_frac020_1000cells/chr*/scent_candidates_chr*.csv`:

| | |
|---|---|
| Candidate pairs | **117,811** |
| Genes | **9,891** |
| Peaks | **46,936** |
| Overlap with the 15,806 LinkPeaks candidates | 5,768 pairs = **36.5%** (42.2% on shared genes) |

So this is not a subset of the LinkPeaks universe. It is 7.5× larger, covers 2.3× more genes,
and disagrees with LinkPeaks about most of the candidate space. **Do not delete these files**
(`TODO.md` §4.3).

Two notes. First, a distance-filtered join over expressed genes × accessible peaks is not a weak
candidate step — it is what LinkPeaks, ArchR, Cicero and SCENT all do. Candidate generation is
not where methods differentiate. Second, the `min_pair_frac` and window values used
(0.02, 100 kb) were chosen for SCENT tractability, not for standalone scoring, and should be
revisited.

### The window-symmetry trap

The committed candidate set is ±100 kb. The SCENT validator is also ±100 kb. Scoring this set
and validating against SCENT would give candidates and validator **identical support**, making
the proximity confound worse than in the current benchmark, with no proximal-removal escape
route — removing links above 50 kb would leave almost nothing testable.

Break the symmetry deliberately. Either generate candidates at 500 kb to match the LinkPeaks
window and treat the 100–500 kb band as explicitly unvalidatable, or accept that a ±100 kb
experiment answers only a proximal question and say so. Do not run the symmetric version and
interpret the result.

---

## 4. Scoring

Keep the score unchanged for v0. Changing candidate generation and scoring simultaneously makes
the result uninterpretable.

$$
S_{pg} = A_{pg} \cdot f_{\lambda}(D_{pg}) \cdot (1 + \alpha T_p)
$$

with \(A_{pg}\), \(D_{pg}\) and \(T_p\) exactly as in `docs/method_report.md` §6–8.

Distance variants to carry forward:

- \(f^{\mathrm{mod}}_{\lambda}(D) = 1 + \lambda(D - 0.5)\) at \(\lambda \in \{0.1, 0.2\}\) — the
  cleaner form, symmetric about \(d_0\)
- \(f^{\mathrm{orig}}_{\lambda}(D) = (1-\lambda) + \lambda D\) at \(\lambda = 0.1\) — for
  continuity with the current results
- A **promoter-dampened or hump-shaped prior** — the one genuinely new variant worth testing,
  since the current family is monotone in proximity and therefore cannot avoid rewarding
  promoter collapse. Something that peaks at intermediate distance rather than at \(d = 0\).

Do not add further \(\lambda\) or \(\alpha\) values without a better validator. The current
sweep already shows the parameter surface is flat where it matters.

Two prerequisites before anything cell-type-specific (§7):

- **Cell-type annotation.** The pipeline clusters at `cluster_resolution: 0.5` and stops. There
  is no marker-based or reference-based annotation step. "Cell-type-specific" currently means
  "unlabelled-cluster-specific", which is not usable output.
- **Metacell or pseudobulk aggregation.** \(A_{pg}\) is a mean over cells of clipped z-score
  products. Splitting the population into cell types reduces the cells per estimate against
  already near-binary ATAC, so cell-type coactivity will be noisier than the global version —
  and the global version is what currently fails to beat `distance_only`. Aggregation is likely
  a precondition, not polish.

---

### Gene- and cell-type-aware TF/motif support

The current TF/motif term should not be carried forward unchanged as a central signal. In this
benchmark it is a peak-level regulatory-potential score:

    peak has motifs for expressed TFs -> peak gets higher score

That construction gives every gene paired to the same peak the same TF/motif score. It cannot
say whether a motif-bearing peak plausibly regulates a particular target gene, and it cannot
distinguish TF programs active in different cell states. In standalone work, the term should be
redesigned from a peak-level score,

\[
T_p
\]

to a peak-gene-cell-state score,

\[
T_{pgc}
\]

where \(p\) is the peak, \(g\) is the candidate target gene and \(c\) is a cell type, metacell
state or latent cell state.

The intended question should change from:

    Does this peak contain motifs for expressed TFs?

to:

    Does this peak contain motifs for TFs that are active in this cell state and plausibly regulate this gene?

A useful standalone TF/motif score should combine three evidence layers:

\[
T_{pgc}
=
\mathrm{motif}_{pc}
\times
\mathrm{TFactivity}_{tc}
\times
\mathrm{TFtarget}_{tgc}
\]

where the components are interpreted as:

| Component | Meaning |
|---|---|
| \(\mathrm{motif}_{pc}\) | The peak contains a motif, motif module or motif grammar active/accessible in cell state \(c\). |
| \(\mathrm{TFactivity}_{tc}\) | The matching TF is active in the same cell state, preferably by motif activity or regulon activity rather than raw RNA alone. |
| \(\mathrm{TFtarget}_{tgc}\) | The target gene behaves like a plausible target of that TF in the same cell state. |

Raw TF RNA expression should not be the main activity proxy. Better activity estimates include
chromVAR motif deviation scores, motif accessibility activity, SCENIC/pySCENIC-style regulon
activity, DoRothEA/VIPER-like TF activity and, if depth allows, footprinting.

The operational question is whether a motif-bearing peak is accessible when the matching TF is
active and the candidate target gene is expressed in the same cell state or metacell context.

This design deliberately separates two things that are collapsed in the current benchmark:

| Feature | Question |
|---|---|
| Peak regulatory potential | Does the peak look like a regulatory element? |
| Gene-specific TF support | Does the motif/TF evidence point to this specific target gene? |

For PBMC-like data, single motifs are likely too noisy. The standalone score should consider
motif families or modules, motif density, motif strength, co-occurring motifs and immune-cell
TF programs such as AP-1, ETS, IRF, NF-\(\kappa\)B, CEBP, RUNX, GATA and TCF/LEF. Ubiquitous
motifs should be downweighted, for example by IDF-like weighting across peaks, collapsing
redundant TF-family motifs, or removing low-information motifs.

External TF-target priors can be used as weak priors rather than hard truth. Candidate sources
include JASPAR motif families, ENCODE/ChIP-derived support where cell-type-relevant, ChEA,
DoRothEA, TRRUST, SCENIC regulons, Perturb-seq and CRISPRi/eQTL evidence where available. The
goal is not to infer a full TF-to-site-to-gene circuit. The goal is to make motif support useful
for peak-gene prioritisation.

The evaluation must include direct motif ablations:

    coactivity only
    coactivity + motif
    coactivity + distance
    coactivity + motif + distance
    motif only
    distance only
    motif + distance only

The key test is stricter than asking whether the full model improves:

    Does motif support improve ranking at fixed distance and fixed coactivity?

If the answer is no, the motif term should remain optional or diagnostic rather than central.

Best-case standalone interpretation:

> A motif helps only when the matching TF is active in the same cell state, the peak is
> accessible in that cell state and the candidate target gene behaves like a target of that TF.

## 5. Baselines and comparison

Same fixed-universe discipline: every method ranked on the identical candidate set.

- `linkpeaks` — run LinkPeaks on the standalone candidate set where possible, so the baseline is
  not advantaged by defining the universe
- `coactivity` alone
- `distance_only` — **mandatory**; the benchmark exists to detect this
- `distance_mod_only`
- `coactivity_tf`, `coactivity_distance`
- `full` at the \(\lambda\) values retained above
- Optional external: ArchR Peak2GeneLinks, Cicero. ArchR was attempted during this project and
  did not produce stable benchmarkable output (`docs/lab_notebook.md`); treat as optional.

### SCENT comparison

SCENT output already exists for 22 autosomes: 52,482 tested rows, 4,758 supporting. If the
standalone candidate set is the same ±100 kb set SCENT was run on, no new SCENT compute is
needed — but see the window-symmetry trap in §3. If candidates move to 500 kb, SCENT coverage
remains 100 kb and the gap must be reported as untested rather than negative.

Reuse the existing support rule (`pvalue_positive`, p ≤ 0.05) and reciprocal-overlap matching
(0.5) so numbers stay comparable to `docs/results_report.md`.

### Distance and proximal controls — carry forward unchanged

1. `distance_only` as an explicit ranking mode
2. Distance-matched enrichment within bins, with bins beyond the validator's window marked
   untested rather than reported as zero
3. Proximal removal at 10 / 25 / 50 kb, with the surviving median distance reported alongside
   every support fraction
4. **Fix the two known defects:** guard the degenerate empty-cell odds ratios so they emit `NA`
   rather than 8.906, and add the missing plot
   (`docs/results_report.md` §6)

---

## 6. Go / no-go criteria

Time-box to two weeks. Write the stop rule down before starting.

**Go** if all of:

1. The standalone candidate set recovers a reasonable share of LinkPeaks' SCENT-supported pairs
   — establishing the candidate step is not silently dropping real links.
2. Within distance bins up to the validator's window, the full model's top-decile SCENT
   enrichment exceeds both LinkPeaks and `distance_only`, at magnitudes comparable to the
   current 5.06 / 3.15 / 2.64.
3. **Tier separation is measurable**: tier 1 links show a materially higher supported fraction
   than tier 3 on the standalone universe. This is the criterion that matters for §7.
4. After proximal removal at 25 kb, and with candidates restricted to the validator's tested
   window, the full model retains an advantage over `distance_only` **and over `coactivity`
   alone**. The current benchmark already clears the `distance_only` bar once the window is
   applied — by +0.28 to +0.32 at \(\delta\) = 25 kb — so that comparison no longer
   discriminates. Beating unadorned coactivity is the bar that does: on this dataset the
   composite score does not beat it in the outermost testable bin.

**No-go** if any of:

1. The full model does not beat `distance_only` within distance bins, or does not beat
   `coactivity` alone.
2. Tiers do not separate on the standalone universe. If tiers do not separate globally they will
   not separate after splitting the cells by type.
3. Candidate generation drops most SCENT-supported pairs.
4. Improvement appears only at parameter settings that increase promoter fraction.

On no-go: stop, and write it up as a second negative result. That outcome still strengthens the
public repository, because it closes the question the current benchmark leaves open — the
asymmetry in value between the two outcomes is what makes the two weeks worth spending.

---

## 7. The actual objective

Worth stating plainly, because it changes what has to be proven.

A twelfth peak–gene linker competing on ranking accuracy against Signac, ArchR, Cicero, SCENT,
SCARlink, CREMA, TRIPOD, SCENIC+, Pando, LINGER and FigR is not a promising target. The score is
the field's consensus ingredient list, and the current evidence does not show it beating a
distance prior.

The more defensible objective is a **cell-type-specific, tiered, calibrated candidate list for
experimental follow-up**. That changes the required claim from

> this ranks better than existing methods

which this evidence cannot support, to

> tier 1 links have measurably higher orthogonal support than tier 3 links, and here is the rate

which is a calibration claim, testable with the harness that already exists, and does not
require beating anyone. Existing tools output score tables; few output a tier with a stated hit
rate and an explicit proximity warning.

Three things must be true for that to be worth building, and only the first is currently known:

1. Tiers separate on a candidate universe the method chose itself — the §6 gate.
2. Cell-type-stratified coactivity is estimable at the available cell numbers — needs the
   metacell work in §4.
3. The output format is something an experimentalist would act on — **not answerable from inside
   the repository.** Before building this, ask one person who runs enhancer perturbations what
   would make them pick up such a list. Perturbation experiments test tens of links, not
   hundreds, so the top of the list must be right — and the top of the list is exactly where the
   proximity confound bites hardest.

## 8. Deliberately out of scope for v0

Deep learning or graph neural network scoring; Hi-C or other chromatin-contact integration; TF
footprinting; multi-tissue or multi-dataset benchmark expansion; causal inference; TF→site→gene
circuit inference (that is CREMA's question — see `docs/competitor_positioning.md`).

One longer-term item does belong on the record, though not in v0: the deeper constraint is that
every validator available here is correlational and derived from the same two matrices. An
orthogonal, distance-decoupled truth set — CRISPRi enhancer-perturbation effect sizes, or
fine-mapped eQTLs — is what would convert any result, positive or negative, into something
defensible. That likely means working in a cell type where such data exists, which PBMC is not.
That is a strategic decision about the evaluation axis, and it should be made before more
scoring code is written.

A future standalone model should be treated as an HPC-scale experiment rather than a local
reranking benchmark. In particular, de novo candidate generation at 500 kb or larger cis
windows, repeated SCENT-style validation, cell-type-stratified analyses and multi-dataset
checks would likely require batch scheduling and larger memory/CPU budgets. The present
release intentionally stops short of that scope.