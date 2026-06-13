# Notes

## for the cli add a module (make linkpeaks) saves to data/results to its own folder where reads it

## Benchmark

## Different distance settings

## Boost?

3. Conditional boosting (best long-term)

Only boost when multiple signals agree:

boost = 1 + β * (distance_score * tf_score)

or

boost = 1 + β * (mul_weigh * tf_score)

Effect:

proximal + TF-supported → boosted
distal but TF-supported → still boosted
proximal but no TF → not overly boosted

👉 This avoids “distance-only domination”

What your current plot already tells you

From your scatter:

vertical spread already exists
that spread is likely dominated by distance

So if you now add boosting on top of distance:

👉 you amplify the same signal twice

What you actually want

Your model should behave like:

Case	Desired behavior
proximal + strong co-activity + TF	strong boost
proximal only	mild boost
distal + strong TF + co-activity	retain or boost slightly
distal + weak signal	penalize

Right now you're good at:

penalizing weak links
promoting proximal plausible ones

You are weaker at:

rescuing strong distal regulatory links
Practical recommendation
Step 1 (do this first)

Keep current model.

Tune:

distance_d0

Goal:

reduce overcompression toward very short distances
Step 2 (then test boosting)

Introduce non-distance boost first:

final_v7 = mul_weigh *
           ((1 - λ) + λ * distance_score) *
           (1 + α * tf_score)

Then optionally:

final_v8 = final_v7 * (1 + β * mul_weigh)
Step 3 (only then test distance boosting)

If you want:

distance_modifier = 1 + λ * (distance_score - 0.5)

Use:

smaller λ (e.g. 0.15–0.2)
and compare carefully
What to watch for immediately

After boosting, check:

Red flags
median top50 distance drops even further (<5 kb)
distal fraction collapses (<10%)
top genes become mostly promoter-driven
ORA becomes narrower (less diverse)
Good signs
distal links with strong TF support reappear
unique genes stay high or increase
ORA stays rich
top-N overlap doesn’t collapse
Bottom line
Yes, naive boosting will mostly boost proximal links
That’s because distance is already your dominant signal

## Note 260418

ArchR was attempted but excluded from the initial benchmark due to environment-specific LSI failures (ArchR 1.0.3, Linux/conda).
### ArchR default IterativeLSI settings were unstable in this environment
(subscript out-of-bounds during LSI). Use a conservative 1-iteration setup
for benchmark robustness; revisit defaults later in containerized pipeline work.

Will be reintroduced in a controlled containerized setup later

## ArchR (Peak2GeneLinks) – benchmarking integration notes

### Status

Attempted integration of ArchR Peak2GeneLinks into the benchmark pipeline.
Pipeline executes end-to-end after multiple adjustments, but output is currently not robust (empty link set under tested configuration).

---

### What works

* Arrow file creation (`createArrowFiles`)
* RNA integration (`addGeneExpressionMatrix`)
* Dimensionality reduction (`addIterativeLSI`, constrained settings)
* Clustering (requires subsampling to avoid pathological runtime)
* Peak calling:

  * `addGroupCoverages`
  * `addReproduciblePeakSet`
  * produces ~62k union peaks
* Peak matrix construction (`addPeakMatrix`)
* Peak2Gene computation step completes without runtime errors

---

### Required adjustments

* Reduced LSI iterations to 1
* Subsampled clustering (`sampleCells`)
* Removed unsupported arguments (`resolution`, `force`) from `addPeak2GeneLinks`
* Added explicit peak workflow (coverage → peaks → peak matrix)
* Relaxed retrieval thresholds:

  * `corCutOff = -1`
  * `FDRCutOff = 1`
  * `varCutOffATAC = 0`
  * `varCutOffRNA = 0`

---

### Current issue

* `getPeak2GeneLinks()` returns **zero links**

  * both with `returnLoops = FALSE` (index-based output)
  * and under relaxed thresholds
* Indicates either:

  * correlations are weak / filtered internally
  * or this configuration does not produce usable links for this dataset

---

### Interpretation

* The ArchR pipeline is operational but **not producing benchmarkable output**
* Issue is not:

  * peak calling
  * matrix construction
* Issue is specifically:

  * **link recovery / usable signal generation**

---

### Practical implications

* ArchR requires:

  * tight coupling between steps
  * version-sensitive behavior
  * non-trivial parameter tuning
* Compared to other methods in this benchmark:

  * higher setup complexity
  * lower robustness in ad hoc execution

---

### Decision (prototype phase)

* ArchR is **deferred from the benchmark comparison**
* Reason:

  * output not reliably obtainable under current lightweight setup
* Will revisit in:

  * containerized environment
  * fully controlled pipeline (versions + parameters fixed)

---

### Next step (optional check)

* Test `returnLoops = TRUE` with relaxed thresholds
* If still empty → stop further debugging

---

### Summary

ArchR integration is partially successful at the pipeline level but currently fails to produce usable peak–gene links for benchmarking. Deferred for later, more controlled implementation.

## ArchR Peak2GeneLinks benchmarking attempt — current outcome

### Summary

ArchR integration was implemented and debugged to the point where the pipeline completed all major preprocessing and peak-to-gene correlation steps. However, under the current lightweight benchmark configuration, ArchR did not yield usable peak–gene links for downstream comparison.

### What worked

The following ArchR steps completed successfully:

* Arrow file creation
* RNA import and addition of `GeneExpressionMatrix`
* IterativeLSI (after switching to a conservative 1-iteration configuration)
* clustering
* group coverages
* reproducible peak set generation
* `PeakMatrix` construction
* `addPeak2GeneLinks`

This means the pipeline itself is operational.

### What had to be changed

Several non-default adjustments were required to make ArchR run at all in the current environment:

* conservative `addIterativeLSI()` settings
* clustering with subsampling
* explicit peak workflow:

  * `addGroupCoverages()`
  * `addReproduciblePeakSet()`
  * `addPeakMatrix()`
* removal of unsupported `addPeak2GeneLinks()` arguments
* relaxed retrieval thresholds in `getPeak2GeneLinks()`

### Current result

Despite successful execution of `addPeak2GeneLinks()`, `getPeak2GeneLinks()` returned zero links even under highly relaxed retrieval settings:

* `corCutOff = -1`
* `FDRCutOff = 1`
* `varCutOffATAC = 0`
* `varCutOffRNA = 0`

This indicates that the issue is not merely downstream filtering or extraction logic. Under the tested configuration, ArchR is not producing a usable benchmark output on this dataset.

### Interpretation

The most likely explanation is that the current ArchR setup is too weak or too coarse to recover stable peak–gene signal in this run. Plausible contributors include:

* conservative LSI/clustering settings
* subsampling
* use of `peakMethod = "Tiles"` rather than a more standard MACS2-based peak set
* dataset/method sensitivity

### Practical decision

For the prototype benchmark, ArchR should be treated as **attempted but excluded** unless one final targeted test is run with:

* `peakMethod = "Macs2"`

If that still returns zero usable links, ArchR should be deferred until later pipeline hardening / containerized implementation.

### Prototype-phase conclusion

ArchR was made operational but did not yield benchmarkable peak–gene links under the current settings. For now, the benchmark should proceed with:

* LinkPeaks
* the reranker
* co-activity baseline
* distance-only baseline
* SCENT

ArchR can be revisited later under a more controlled and better-tuned setup.


## ArchR MACS2 follow-up

A final ArchR retry was attempted using `peakMethod = "Macs2"` to test whether the empty-link result from the tile-based peak workflow was due to coarse peak definition.

### Outcome

This retry did not proceed to peak calling because MACS2 was not available in the environment:

* ArchR searched for MACS2 in `$PATH`
* also checked `pip` / `pip3`
* no executable was found

Resulting error:

> `Could Not Find Macs2!`

### Interpretation

This means the MACS2-based ArchR follow-up was not a biological negative result. It was a missing-dependency failure.

### Practical takeaway

At this point, ArchR has required multiple environment-specific fixes and still has not produced a clean benchmark output under the current setup. The MACS2 retry adds further evidence that ArchR should be deferred until later pipeline hardening / containerization.

### Current ArchR status

* Tile-based ArchR workflow: runs, but returns zero usable peak–gene links
* MACS2-based retry: blocked by missing external dependency
* Conclusion: ArchR excluded from the prototype benchmark for now

“A constrained ArchR prototype pipeline executes under the current environment, but does not yet produce benchmarkable peak–gene links. The result should be treated as an integration status, not a comparative performance result.”

ArchR integration was attempted, but under the current prototype environment it required multiple workflow-specific adjustments and still did not yield robust benchmarkable output. ArchR was therefore excluded from the prototype comparison and deferred to later containerized pipeline work.

## Per cell type scores and make cell type annotations + cell type groups + all + comparissons
### Here I would think I will annotate the cell types and probably make a noise reduction by making metacells (sort of pseudobulks) like in previous repo
### so make the scoring based on metacells? Test if it is "more stablile" or is there any justification for it??? Maybe



#  transcript-derived TSS table

The recent revision improves both the biological validity and the evaluation integrity of the reranking pipeline. 
The most important change is the replacement of the old gene-coordinate collapse logic with a dedicated transcript-derived TSS table. 
Previously, gene TSS positions were approximated by aggregating coordinates from the Signac annotation object, 
which could mix multiple transcript-derived positions and distort peak-to-gene distances. The updated implementation now derives 
one explicit TSS per gene directly from EnsDb transcript annotations, with a preference for protein-coding and longer transcripts. 
Because distance is a central component of the scoring function, this change makes the distance prior substantially more trustworthy 
and reduces the risk that links were promoted or penalized due to an artifact of annotation collapse rather than real genomic proximity.

The second major improvement is the correction of the baseline evaluation framing. Earlier, some “baseline” summaries were effectively
 being computed within the reranked candidate subset, which made the baseline comparison cleaner than it really was. The revised script
 now separates the full LinkPeaks baseline from the reranked candidate table. Full baseline links are deduplicated once, retained as their own object, 
 and used for baseline ORA, overlap summaries, and baseline distance calculations. In addition, a dedicated baseline_dist table is now built so that
 baseline distance summaries are computed from the true full baseline rather than from the restricted reranking results. This makes the comparisons 
 against the reranked model much more honest: the baseline is now a real baseline, not a conditional subset of the method’s candidate space.

Several robustness and reproducibility fixes were also added. The script now sets a random seed from the CLI, which improves reproducibility
 for steps such as clustering and neighborhood construction. Fragment index validation was added, so the workflow now fails early if the required
 .tbi file is missing instead of breaking later in a less interpretable way. The motif-score rescaling step was also made safer by replacing
 dimension-dropping apply() behavior with a more stable column-wise reconstruction. In parallel, some unnecessary or misleading logic was removed or 
 cleaned up, including the pointless save/reload cycle of the Seurat object, the duplicate motif_tf_score field, and old unused TSS code.

The net effect of these changes is that the script is now much better aligned with the scientific questions it is trying to answer. 
The reranking formula itself did not fundamentally change, but the biological inputs to the distance prior are now more defensible, 
and the baseline-vs-reranked comparisons are now methodologically sound. In practice, this means future results from the pipeline should
 be interpreted as more reliable: if the reranked model still outperforms the baseline after these corrections, that improvement will carry substantially
more weight. The script remains a research prototype rather than a production workflow, but it is now markedly stronger as a publication-grade analysis foundation. 


### De novo model

Full Model Plan (Post-Holiday)
Goal

Move from reranking LinkPeaks candidates → standalone peak–gene linking model

Key shift:

control candidate generation
apply own scoring model
reduce dependence on external methods
Core Design
1. Candidate generation (replace LinkPeaks)

For each gene:

same chromosome
within cis window (start: 500 kb)
filter:
gene expression prevalence
peak accessibility prevalence

Output:

candidate pairs: (gene, peak)
2. Cluster-aware structure
cluster cells (coarse is fine initially)
optionally annotate later
treat clusters as independent regulatory contexts
3. Metacell construction (within cluster)

Within each cluster:

group cells → metacells
aggregate:
RNA expression
ATAC accessibility

Purpose:

reduce sparsity
stabilize co-activity signal
improve distal link detection
4. Scoring (per cluster)

For each cluster k:

Core model:

A_pgk = coactivity (metacell-level)
D_pg  = distance prior
T_pgk = TF/motif support

S_pgk = A_pgk × D_pg × (1 + α T_pgk)

Where:

coactivity = positive co-expression/accessibility coupling
distance = smooth decay (not hard cutoff)
TF = motif × TF expression (cluster-aware)
5. Combine across clusters

For each (peak, gene):

Recommended first version:

S_pg = max_k S_pgk
best_cluster = argmax_k S_pgk

Optional later:

number of supporting clusters
weighted average variants
6. Output structure

Per link:

peak
gene
global_score
best_cluster
cluster_score
distance_bp
tf_score
(optional: n_clusters_supported)
Benchmark positioning

Important:

This model operates on its own candidate universe
Not directly comparable by pair overlap unless peak sets are harmonized

Comparison strategy:

distance distribution
distal fraction
gene diversity
enrichment (ORA)
top-ranked behavior
Known issues (from current work)
Peak universes differ across methods → breaks pair overlap
SCENT required fallback to native candidates
ArchR unstable → exclude for now

Conclusion:

Do not enforce shared candidate universe for all methods at this stage
Immediate next steps (after return)

APPLY:

implement cis candidate generation
implement metacells within cluster
port existing scoring formula to metacell level
implement cluster-wise scoring + max aggregation

SKIP:

complex bias correction
peak harmonization across methods
ArchR reintegration
advanced modeling (deep learning, GNN, etc.)
Expected improvements

If working correctly:

more stable scores than raw-cell version
better distal link recovery
clearer cluster-specific regulation
improved biological coherence
One-line summary

Build a standalone peak–gene linking model using cis candidates, cluster-aware metacells, and your current scoring formula applied per cluster, then aggregate across clusters.

That makes sense.

That is a good progression:

Phase 1

Prove the current method has value:

benchmark current reranker
clean up evaluation
show it improves over baseline
Phase 2

Turn it into a more complete pipeline:

start from raw FASTQ
standard preprocessing
candidate generation
reranking
evaluation
Phase 3

Move toward standalone inference:

replace LinkPeaks candidate generation
use your own candidate builder
then later add cell-type-aware / metacell logic

That is a sane order.

Important point

Starting from raw FASTQ does not by itself make it more de novo.

It makes it a more complete and usable pipeline, which is good for:

repo quality
reproducibility
publication packaging

But the real step toward de novo is still:

own candidate generation instead of LinkPeaks top-K

So I would separate those mentally:

FASTQ pipeline = engineering/completeness upgrade
own candidate generation = methodological upgrade

Both matter, but they are different.

Best practical route

If benchmarking looks good, I would do:

1. Repo / pipeline upgrade
raw FASTQ → counts/fragments
standardized run config
reproducible outputs
2. Candidate generator v2

Replace LinkPeaks candidates with:

same chromosome
within distance window
minimum gene expression prevalence
minimum peak accessibility prevalence

That is the easiest first standalone version.

3. Keep LinkPeaks mode too

Do not throw it away.

Support both:

candidate_mode = "linkpeaks"
candidate_mode = "window"

That gives you:

backward compatibility
cleaner benchmarking
a path to show how much your own candidate generation changes results
Why this is a good journey

Because each stage has value on its own:

current reranker = useful method prototype
raw FASTQ pipeline = strong portfolio/repo
own candidate generation = more novel method
cell-type/metacells = stronger biological specificity

That is a real roadmap, not random feature creep.

One caution

Do not jump from current reranker straight into:

raw FASTQ
own candidate generation
cell-type specificity
NB model

all at once.

That would make it impossible to know what improved anything.

Better sequence
finish current benchmark
freeze reranker version
build end-to-end pipeline
swap candidate generator
benchmark again
only then move to cell-type-aware extension
Bottom line

Yes — that is a good path.

If current benchmarking is reasonable, an end-to-end pipeline from raw FASTQ is a strong next step. Just keep in mind that:

FASTQ support makes it more complete
own candidate generation makes it more de novo

Those are the two separate upgrades.

 

The good sign in your case

A lot of weak methods fail the smell test early. Yours does not.

Your outputs already look:

interpretable
biologically plausible
tunable
benchmarkable

That is exactly the profile of something worth continuing.

Bottom line

Yes — this could absolutely become a reasonably good method.

Not because it is flashy, but because it is:

scoped correctly
evidence-integrating
interpretable
already showing useful behavior

The next thing that matters is not more intuition. It is whether it beats:

baseline
distance-only
simplified ablations

If it does, then you likely have something real.


# 260417

Thought for a couple of seconds

Yes — the screenshot matches the code’s intent.

Your current scoring function is essentially:

S
pg
	​

=(
∣C∣
1
	​

c∈C
∑
	​

max(z
gc
RNA
	​

,0)max(z
pc
ATAC
	​

,0))⋅[(1−λ)+λ
1+(d
pg
	​

/d
0
	​

)
2
1
	​

]⋅(1+αt
p
	​

)

So the pipeline is already in the right shape for a prototype reranker:

term 1 = joint RNA/ATAC support
term 2 = distance prior
term 3 = TF/motif support

That is a clean starting point for a gradual move toward a de novo model.

What I would do next

Do not jump straight to a new model class yet. First, make this reranker into a proper benchmark harness.

Phase 1: lock down the current reranker

Before changing the model, test whether each term is actually helping.

Run these ablations:

coactivity only
S=mul_weigh
coactivity + distance
S=final_v5
coactivity + TF
S=mul_weigh(1+αt
p
	​

)
distance only
TF only
full model
S=final_v6

This tells you whether the distance term and TF term add signal or just add plausibility bias.

Benchmarking targets to add

You need at least four benchmark layers.

1. Internal ranking stability

Test whether results are stable across:

random seeds
downsampling cells
downsampling fragments
different clustering resolutions
different candidate_top_k values

Useful metrics:

top-100 / top-500 overlap
Spearman rank correlation
stability of promoted links
stability of top genes per cell type

If the ranking changes wildly, the model is not ready for de novo expansion.

2. Biological plausibility

Your current ORA is fine as a first pass, but it is weak as a primary benchmark.

Add:

promoter vs distal composition
known marker-gene proximal peaks
motif consistency by lineage
TF–target coherence

Examples:

in hematopoietic data, GATA / SPI1 / CEBP family peaks should enrich near expected lineage genes
in neural data, NEUROD / SOX / OLIG families should behave sensibly

This is still soft evidence, but much better than GO alone.

3. External truth sets

This is the main step.

Benchmark against datasets with some orthogonal regulatory evidence:

CRISPR perturbation enhancer–gene maps
promoter capture Hi-C / pcHi-C
ABC or EpiMap enhancer-gene references
validated enhancer databases
eQTL-linked regulatory regions, where applicable

Do not treat any of these as absolute truth. Treat them as noisy positives and compare enrichment.

Recommended metrics:

AUROC
AUPRC
recall@K
enrichment of supported links in top-N vs baseline
odds ratio in top decile vs all candidates

For sparse truth sets, precision-recall and recall@K matter more than AUROC.

4. Baseline comparisons

Your real baseline set should be larger than:

LinkPeaks
distance only

Also compare against:

nearest-gene
correlation only
cicero / co-accessibility-style baseline if possible
ABC-style score approximation if you can construct one
simple logistic or rank-sum ensemble using the same features

If a linear/logistic baseline beats your handcrafted score, that is useful information.

Best next model step

APPLY: move from hand-tuned multiplicative scoring to a supervised or weakly supervised ranker.

Good next options:

Option A: logistic regression / gradient boosting on link features

Build one row per candidate link with features like:

LinkPeaks score
RNA–ATAC coactivity metrics
distance
GC content
peak accessibility mean
gene expression mean
motif score
TF score
promoter indicator
same-cluster specificity
pseudobulk correlation
cell-type-specific correlation

Then train:

logistic regression
XGBoost / LightGBM ranker
pairwise ranking model

This is the fastest serious upgrade.

Option B: pseudobulk-aware de novo score

Instead of cell-level z-product only, compute within-cluster or pseudobulk summaries:

peak accessibility by cluster
gene expression by cluster
correlation across clusters
specificity matching

That often works better than cell-level noise for enhancer–gene ranking.

Option C: graph-based model later

Once features and benchmarks are solid, then move to:

bipartite peak–gene graph
GNN or message passing
latent variable model over peak/gene/cell-type nodes

SKIP this for now unless you already have strong training labels and reproducible feature wins.

What to fix in the current prototype before scaling
1. Separate candidate generation from ranking

Right now LinkPeaks is the gateway.

For a de novo model, create candidate sets independently, such as:

all peaks within 1 Mb of TSS
all peaks within TAD if available
promoter and distal handled separately

Then let the model rank them.

This is the single most important architectural shift.

2. Add cell-type-aware scoring

Your current TF score is global. Better:

compute TF expression per cluster
compute motif support per cluster
score peak–gene links within matched cellular contexts

That will help a lot.

3. Stop using one transcript/TSS per gene as the long-term design

For prototyping it is fine. For a better model:

use canonical transcript set or gene-level promoter windows
maybe multiple TSSs per gene
distinguish promoter-proximal vs distal links
4. Add negative controls

You need deliberate false links:

cross-chromosome links
far-distance matched negatives
shuffled peak–gene pairs preserving marginals

These are critical for debugging whether the score is learning biology or just distance.

Concrete benchmark plan

I would structure the next work like this.

Sprint 1
freeze current pipeline
add ablation runs
add stability metrics
add candidate-independent nearest-gene baseline
add shuffled negatives
Sprint 2
build feature table for every candidate link
run logistic regression and XGBoost baselines
compare against current formula with recall@K and PR-AUC
Sprint 3
add external truth sets
benchmark by distance bins:
promoter
0–50 kb
50–200 kb
200–500 kb
500 kb+

This matters because many models look good only because they win near promoters.

Sprint 4
replace LinkPeaks candidate generation
generate full cis windows
retrain/retest
only then call it de novo-ish
What success should look like

A strong next milestone is not “new model.”
It is:

full model beats LinkPeaks and distance-only on recall@K
improvements persist after ablations
gains hold across seeds and downsampling
gains are strongest in distal bins, not only promoter bins
a simple learned ranker on your features performs at least as well as the hand-designed formula

If you get that, you have a real platform to build a de novo model from.

Recommendation

APPLY:

turn this into a benchmarked reranking framework first
add ablations, negatives, external truth sets, and feature-based baselines
then move to candidate-independent ranking

SKIP:

GNNs, VAEs, or complex de novo architectures right now
claiming de novo inference while LinkPeaks still defines the candidate universe

Your current prototype is good enough to justify the next step. The immediate priority is not a fancier score. It is a harder benchmark.

I can help you turn this into a concrete benchmark matrix with exact metrics, baselines, and output tables.


I know. I am doing this shit you list

Good — then you’re on the right track. Next step is tightening execution so you don’t waste cycles.

What to focus on now (high ROI)
1. Lock a minimal benchmark suite

Do not keep expanding. Fix a small, repeatable set:

Datasets (2–3 max)

1 clean multiome (good depth)
1 noisier dataset
optionally 1 with known regulatory signal (e.g., hematopoietic)

Truth sources (pick 1–2)

ABC / EpiMap-style enhancer–gene
pcHi-C or CRISPR if available

Metrics (don’t exceed these)

recall@50 / @100
AUPRC
enrichment vs distance baseline (odds ratio top decile)

Anything beyond this is noise at this stage.

2. Make your current model diagnosable

Right now your score is clean but opaque.

Add these columns and inspect distributions:

coactivity = mul_weigh
distance_score
tf_score
final_v6

Then compute:

correlation matrix of features
feature contribution to top-N (mean feature value in top 100 vs all)

If TF or distance barely shifts ranking → drop or rescale.

3. Run controlled ablations once, not repeatedly

Run exactly this matrix and freeze it:

model	components
M1	coactivity
M2	coactivity + distance
M3	coactivity + TF
M4	full
M5	distance only
M6	LinkPeaks

Export one CSV summarizing:

recall@K
AUPRC
distal recall (>50kb)

Do not re-run unless you change features.

4. Add one simple learned baseline

This is critical.

Build a table:

features = [
  link_score,
  mul_weigh,
  distance_score,
  tf_score,
  motif_score
]

Train:

logistic regression
or XGBoost (depth 3–5 max)

If this beats your hand formula → switch immediately.

APPLY if better
SKIP hand-tuning further

5. Check distance leakage

Very common failure mode.

Do this:

bucket links by distance:
0–10kb
10–50kb
50–200kb
200–500kb

500kb

Then compute recall@K per bucket.

If all gains are in <50kb → your model is just a distance proxy.

6. Kill weak components early

Use this rule:

If adding TF score improves <5% recall@100 → drop or rethink
If distance term dominates ranking (>0.8 correlation with final score) → reduce λ
If LinkPeaks still dominates → your reranker isn’t adding signal
7. Prepare for de novo (but don’t switch yet)

You only move when:

learned baseline ≥ your current model
gains persist in distal bins
results stable across seeds

Then:

Replace candidate set

all peaks within 1 Mb of TSS
optionally TAD-restricted

Everything else stays the same initially.

Quick sanity checklist

You’re ready to move forward if:

 full model beats LinkPeaks on recall@100
 beats distance-only clearly
 gains not confined to promoters
 feature importance is not dominated by one term
 results stable under downsampling

If any of these fail → fix before expanding.

Bottom line

You don’t need more ideas right now. You need:

a fixed benchmark
one learned baseline
strict ablation discipline

Do that once, cleanly, and you’ll know exactly whether your formulation is worth scaling into a de novo model or not.

If you want, next step I can:

sketch the exact feature table + training code (R or Python)
or define a minimal evaluation script so you stop rewriting analysis each run

One clean sentence for your writeup

A distance-only ranking collapsed to promoter-proximal links (median ~11 bp) and yielded no significant enrichment, 
indicating that genomic proximity alone is insufficient to recover biologically coherent regulatory programs.

# Benchmarking Notes

What to actually run (practical set)
Core (DO THESE)

APPLY

LinkPeaks
ArchR Peak2GeneLinks
SCENT

These are enough to establish credibility.

Add 1 strong modern method

Conditional APPLY

scMultiMap → only if setup time < ~1–2 days

Reason:

It’s strong, but setup/integration cost is non-trivial
If it slows you down, it’s not worth it right now
Simple baselines (MANDATORY)

APPLY

distance-only ✅ (you already have it)
nearest-gene
correlation-only (your mul or mul_weigh ≈ already this)

These are critical for reviewers.

Optional / secondary

APPLY (optional)

Cicero → include, but label clearly as co-accessibility

If time is tight:
SKIP Cicero initially

Clean benchmark tiers
Minimum viable (what you should do now)
distance-only
nearest-gene
LinkPeaks
ArchR
SCENT
your model

This is already a publishable benchmark set.

Expanded (only if time allows)
scMultiMap
Cicero
What this person got right

Correct:

distinguishing same-task vs related-task (Cicero)
insisting on simple baselines
prioritizing multimodal methods (SCENT, scMultiMap)
What they slightly overreach on
Including everything at once
→ slows you down massively
scMultiMap as “must-have”
→ not true for first pass
Your current position (based on your results)

You already have:

strong improvement over LinkPeaks (ORA ↑, diversity ↑)
distance behaving correctly (not dominating)
TF adding signal (ORA 34 → 52 = real effect)

You are ready for benchmarking now, not after building a massive panel.

Execution plan (do this)
Step 1 (immediate)

Run:

nearest-gene
correlation-only (reuse mul_weigh)
LinkPeaks (done)
your model
Step 2

Add:

ArchR Peak2GeneLinks
SCENT
Step 3 (optional)

Try:

scMultiMap
Bottom line

APPLY:

LinkPeaks
ArchR
SCENT
distance-only
nearest-gene
correlation baseline

OPTIONAL:

scMultiMap
Cicero

SKIP (for now):

trying to include everything at once

If you want, I can give you a minimal benchmarking table schema + plotting plan so you don’t waste time structuring results.

Step 2

Add:

ArchR Peak2GeneLinks
SCENT

Ok. I will do these for now. And if those are somehow reasonable try with breastcancer and put scMultiMap in backlog if it is heavy

Good plan. This is the right level of scope.

What you should do next (tight)
1. Implement ArchR + SCENT

Run them on the exact same dataset + gene universe.

Do NOT change:

filtering
genome build
gene naming

Keep everything aligned.

2. Standardize outputs immediately

Before you go too far, enforce this structure for every method:

peak | gene | score | method

You will save a lot of pain later.

3. Evaluate with the same metrics

For each method compute:

top-N genes (same N as yours)
ORA
median distance (top50)
distal fraction (>50kb)
overlap with your model (optional)

No custom metrics per method.

4. What you are looking for

You want to answer:

A. Do you beat LinkPeaks?

(already yes)

B. Do you match or beat ArchR / SCENT?

Outcomes:

You win clearly → strong result
You match them → still good (simpler method)
You lose slightly but add interpretability → still defensible
You lose badly → revisit model (don’t move to new dataset)
5. Only after that → new dataset

Your plan is correct:

“if those are reasonable → try breast cancer”

Do NOT jump earlier.

6. scMultiMap handling

Your idea is correct:

put in backlog
only implement if:
your method is already competitive
and you need a stronger comparator
One practical warning

ArchR and SCENT can differ in output format:

ArchR → often returns peak–gene links already aggregated
SCENT → may output statistical associations per pair

You must convert both into:

peak | gene | score

before ranking.

Minimal success criterion

After Step 2, you should be able to say:

“Our method improves over LinkPeaks and performs comparably or better than ArchR/SCENT in enrichment while maintaining biologically plausible distances.”

If you can say that, you’re in a strong position.

Bottom line

Proceed exactly as you described.

ArchR + SCENT now
Breast cancer only after validation
scMultiMap later if needed

This is efficient and defensible.

# Distance lambda

1. Key result (most important)
Overlap: λ = 0.3 vs λ = 0.5
Top100 pair overlap: 90%
Top100 gene overlap: 84%
Top200 pair overlap: 182 / 200 (~91%)

👉 Translation:

λ = 0.5 is mostly a reordering of the same links, not discovering new biology.

2. What λ is actually doing (now proven)
λ = 0.5
pushes already good links up
removes distal ones from top ranks
increases enrichment by tightening ranking
λ = 0.3
keeps same core links
preserves distal ones in top ranks

👉 So:

This is ranking refinement, not signal discovery.

3. Why ORA increases at λ = 0.5

Because:

you concentrate on:
immune genes
proximal regulatory structure
you remove:
noisy distal links

So enrichment improves.

But:

It is partly selection bias toward easier biology (proximal)

4. The decisive argument (this is your paper logic)

You can now say:

“Increasing distance regularization (λ = 0.5) increases enrichment (65 vs 52 ORA terms), but does not substantially change the selected links (90% overlap in top100), instead shifting ranking toward promoter-proximal interactions (distal fraction 2% vs 10%).”

That is strong.

5. Final decision (no ambiguity)
KEEP as main model

👉 λ = 0.3

Because:

same core biology as λ = 0.5
preserves distal links (54 vs 9)
avoids promoter collapse
still strong enrichment
REPORT as sensitivity

👉 λ = 0.5

Position:

improves enrichment
shifts toward proximal regulation
does not change link identity substantially
REJECT

👉 λ = 0.1

6. One more important insight (you now proved this)

From overlap + metrics:

Distance is primarily controlling ranking, not discovery

This is actually good:

your model finds links via signal (mul + TF)
distance refines ordering

That’s exactly how it should behave.

7. You are done with tuning

Do NOT:

try more λ
try more α

You already mapped the behavior.

8. Next step (clear)
APPLY now
LinkPeaks
ArchR Peak2GeneLinks
SCENT
After that
second dataset (breast cancer)
Bottom line
λ = 0.5 looks “better” only because it compresses ranking
λ = 0.3 gives the same biology with better regulatory balance
your model is stable across λ → very strong result

Move on to benchmarking.


# alpha TF

What your results actually say
α = 0.1 → weak TF
ORA: 40 (bad)
correlation: highest (0.733) → basically baseline-like
TF influence too weak → not useful

➡️ SKIP

α = 0.3 → balanced
ORA: 52
good distance behavior
high overlap with baseline (~94/100)

➡️ Solid baseline

α = 0.5 → strongest signal
ORA: 56 (best)
gene diversity highest (75)
TF signal clearly stronger:
median_tf_top100 = 0.255 (highest)
high-dist TF also highest
BUT:
correlation drops (0.721)
slight shift in rankings (overlap 91/100 vs baseline)

➡️ BEST overall

Key observation (this matters)

From your overlap:

0.3 vs 0.5:
91/100 pairs overlap
~87% gene overlap

This is critical:

α=0.5 is NOT inventing a new solution — it is refining the same solution with better biological prioritization

That’s exactly what you want.

Should you try α > 0.5?
What will happen if you increase α further

Predictable behavior:

TF dominates ranking
you start selecting:
motif-heavy peaks
possibly noisier / less supported links
correlation drops further
ORA may:
plateau
or degrade
Decision
APPLY
α = 0.5 → use as main model
OPTIONAL (one test only)
α = 0.7
only to confirm saturation
expect:
ORA ≈ same or worse than 0.5
correlation ↓ further
SKIP
anything > 0.7
dense grid search (waste of time)
Final model recommendation

Use:

lambda_distance = 0.3
alpha_tf = 0.5

This is now empirically supported by:

best ORA
stable distance profile
minimal disruption of core ranking
improved TF signal where it matters (distal + high-tier)
What to do next (priority order)
1. LOCK parameters

Stop tuning α and λ.

2. Run comparators

Proceed with:

LinkPeaks (baseline)
ArchR Peak2GeneLinks
SCENT
3. Generate figures (your script is ready)

Use:

meaningful distal plots
custom genes
4. Optional sanity check

Look at:

top promoted links (α=0.3 → α=0.5)

If those are biologically sensible → you're done.

Bottom line
You already found the sweet spot.
α = 0.5 is your final setting.
Don’t over-optimize — move to benchmarking.

# Notes 260422

## SCENT benchmark attempt — runtime outcome

A prototype SCENT benchmark run was executed on the native cis candidate set after shared-universe intersection failed because the SCENT peak set and the LinkPeaks/reranker peak set had zero exact peak overlap.

### Candidate scale

* genes passing prevalence: 8475
* peaks passing prevalence: 28894
* cis candidate pairs after distance filter: 195426

### Runtime behavior

The run was computationally active and healthy throughout:

* 4 worker processes remained CPU-saturated
* memory use was high but stable
* no major swap thrashing
* no evidence of process hang

However, total runtime became impractically long (~75 hours) for a prototype benchmark configuration.

### Interpretation

This is not a SCENT installation/runtime failure. It is a scaling/configuration failure: the current all-cells, ~195k-pair setup is too large to serve as the default benchmark run.

### Decision

The run was stopped and SCENT will be restructured later using a more practical execution strategy, likely:

* per-chromosome runs
* smaller cis windows
* stronger candidate prefiltering
* or smaller biological subsets

### Takeaway

SCENT is functional, but the current benchmark configuration is too expensive to be useful.
Future SCENT benchmarking should focus on tractable subproblems rather than one monolithic all-cells run.

## Post-Holiday Restart Plan (Peak–Gene Benchmark)

### Goal

Quickly regain context and resume development without re-debugging old issues.

---

### 1. Re-establish current state (1–2 hours)

* Re-read:

  * benchmark script (especially SCENT + ArchR blocks)
  * latest notes on ArchR exclusion and SCENT runtime
* Confirm:

  * inputs (h5, fragments, outputs)
  * working methods: LinkPeaks, Reranker, Coactivity, Distance
* Do **not** rerun heavy jobs yet

---

### 2. Fix known structural issues (high priority)

#### Peak mismatch (critical)

* Problem:
  SCENT peaks ≠ LinkPeaks/reranker peaks → zero overlap
* Action:

  * implement **GRanges overlap mapping**
  * use **reciprocal overlap ≥ 0.5**
* Do not use exact string matching anymore

#### Benchmark logic

* SCENT:

  * run on **native candidate set**
  * compare via summaries, not exact pair overlap (initially)
* ArchR:

  * keep **excluded (attempted)** unless fully stabilized

---

### 3. Redesign SCENT execution (core fix)

#### Replace monolithic run

SKIP:

* ~200k pairs
* all cells
* single run

APPLY:

* **per chromosome runs**

  * split candidate set by chr
  * run SCENT per chr
  * merge results

Optional:

* smaller cis window (100–200 kb for first pass)
* test on subset before full run

---

### 4. Validate SCENT output (first successful run)

After first successful run:

* check:

  * number of links
  * score distribution
  * runtime
* compare:

  * gene-level overlap with LinkPeaks / reranker
  * distance distribution
  * ORA behavior

Do not focus on exact pair overlap yet

---

### 5. Benchmark framing (clean)

Include:

* LinkPeaks
* Reranker
* Coactivity
* Distance
* SCENT (native candidate space)

Exclude:

* ArchR → **attempted, not included**

---

### 6. Next development direction (after benchmark stabilizes)

* integrate your model:

  * cluster-aware / metacell-based linking
  * add TF/motif features
  * bias correction
* treat current reranker as baseline

---

### Key takeaways to remember

* Peak mismatch ≠ method failure
* SCENT scaling issue = configuration problem
* ArchR issue = integration complexity, not biology
* Reranker = currently strongest stable component

---

### First command to run after holiday

Do this before anything heavy:

```r
length(intersect(unique(cand$peak), unique(candidate_universe$peak)))
```

This confirms whether peak alignment is fixed.

---

### Stop conditions (important)

If:

* runtime > few hours for test runs
* memory unstable
* candidate count explodes

→ reduce problem size immediately

---

### Working principle

Start small → validate → scale
Not the other way around

## Full Model Restart Plan (Own Cis-Window Candidate Model)

### Main next-method step

Move from:

* reranking LinkPeaks candidates

to:

* **own candidate generation + own scoring**

This is the real methodological upgrade.

---

### 1. Define own cis candidate universe

For each gene:

* same chromosome
* peaks within fixed cis window
* start with:

  * **100–500 kb**
* later sensitivity:

  * 1 Mb if needed

Filters:

* minimum gene expression prevalence
* minimum peak accessibility prevalence

Output:

```text
(gene, peak) candidate table
```

This replaces LinkPeaks as the entry gate.

---

### 2. Keep current scoring logic first

Do **not** redesign the formula immediately.

First standalone version should reuse current logic:

```text
score = coactivity × distance × TF
```

Meaning:

* coactivity = main signal
* distance = soft prior
* TF/motif = biological support

Goal:

* prove current scoring still works when candidate generation is no longer inherited from LinkPeaks

---

### 3. First standalone version

Implement:

#### candidate mode

* `candidate_mode = "window"`

#### scoring inputs

* RNA matrix
* ATAC matrix
* transcript-derived TSS table
* motif / TF support

#### output

* ranked peak–gene links from your own cis window universe

This is the simplest real full-model milestone.

---

### 4. Then add cluster-aware version

After basic cis-window model works:

* cluster cells
* annotate later if needed
* create metacells within cluster
* compute coactivity on metacells instead of raw cells

Within cluster k:

```text
S_pgk = A_pgk × D_pg × (1 + α T_pgk)
```

Then aggregate across clusters, likely:

```text
S_pg = max_k S_pgk
```

Store:

* best cluster
* global score
* maybe number of supporting clusters

---

### 5. Correct development order

APPLY:

1. own cis candidate generation
2. run current formula on that candidate set
3. benchmark standalone version
4. only then add metacells / cluster-aware scoring

SKIP:

* changing candidate generation and scoring logic at the same time
* adding advanced bias correction before basic standalone version works

---

### 6. What to check first once implemented

* candidate count per gene
* promoter vs distal balance
* score distribution
* whether top links collapse to nearest promoter
* whether ORA / gene diversity stay strong

Main question:

> does current scoring still work when freed from LinkPeaks candidate space?

---

### 7. Immediate post-holiday coding target

Implement a clean standalone pipeline with:

* transcript-derived TSS table
* cis-window candidate builder
* prevalence filters
* current formula
* ranked output

That is the real next milestone, more important than further benchmark plumbing.

---

### One-line reminder

The next real method step is **not** more reranking tweaks — it is replacing LinkPeaks candidate generation with your own cis-window candidate model.


## 260607 SCENT chr22 smoke-test summary

## SCENT chr22 smoke-test summary

SCENT is functional on this dataset. The earlier long runtime was mainly a scale/configuration issue, not an installation or candidate-formatting failure.

The key fix was to subset the RNA and ATAC matrices to only the genes and peaks present in the SCENT candidate table before calling `CreateSCENTObj()`.

| Test | Chr | Cells | Window | Min frac | Candidate pairs | Genes | Peaks | Output rows | Runtime |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Smoke test | chr22 | 500 | 25 kb | 0.20 | 82 | 50 | 76 | 82 | completed |
| Scale-up test | chr22 | 1000 | 100 kb | 0.20 | 113 | 53 | 94 | 113 | 16.71 min |

For the 1000-cell / 100 kb run, SCENT returned `gene`, `peak`, `beta`, `se`, `z`, `p`, and `boot_basic_p`. The ranking score was computed as:

`-log10(boot_basic_p) * sign(beta)`

The output had 113 rows, with 59 positive and 54 negative scores.

Next step: integrate the matrix-subsetting fix into `benchmark_methods.r` and run a real chr22 benchmark with LinkPeaks, Reranker, Coactivity, DistanceOnly, and SCENT. Skip full-genome SCENT and ArchR for now.


## 260613

## Benchmark summary: Reranker vs LinkPeaks, Coactivity, DistanceOnly, and SCENT

Current benchmark compares:

| Method       | Role                                      |
| ------------ | ----------------------------------------- |
| LinkPeaks    | baseline peak–gene links                  |
| Coactivity   | RNA/ATAC coactivity-only baseline         |
| DistanceOnly | nearest/proximal control                  |
| SCENT        | external comparator from chromosome sweep |
| Reranker     | current full model                        |

### Main result

The reranker produces the strongest GO BP enrichment among the tested methods. Its top genes are enriched for immune, lymphocyte, T-cell, leukocyte, and adaptive immune response terms, which is biologically plausible for the current dataset.

This is not explained by distance alone. DistanceOnly collapses to promoter-proximal links and produces little/no useful enrichment. The reranker also has low exact-pair overlap with DistanceOnly, so it is not just a nearest-gene ranking.

### Distance behavior

| Method       | Distance behavior                                    |
| ------------ | ---------------------------------------------------- |
| DistanceOnly | promoter/TSS collapse                                |
| Reranker     | strongly proximal, but not identical to DistanceOnly |
| SCENT        | intermediate, mostly within 100 kb                   |
| LinkPeaks    | broadest / most distal-heavy                         |
| Coactivity   | broad / distal-heavy                                 |

The reranker strongly reduces distal links compared with LinkPeaks and Coactivity. This appears to improve biological coherence, but it is also the main caveat: the model may be over-regularizing toward proximal links.

### Overlap behavior

The reranker overlaps substantially with LinkPeaks and Coactivity, but only weakly with DistanceOnly. This suggests the reranker is still using coactivity/baseline signal, while reordering links with distance and TF/motif support.

SCENT overlaps only modestly with the reranker and gives a different ORA profile dominated by translation/ribosome/rRNA terms. Treat SCENT as a useful comparator, not as a gold standard.

### Interpretation

Current result:

**Reranker = LinkPeaks/coactivity-informed ranking with stronger distance and TF/motif regularization.**

This is a promising pilot result because:

| Check                                       | Status           |
| ------------------------------------------- | ---------------- |
| Beats distance-only in biological coherence | yes              |
| Improves ORA over LinkPeaks / Coactivity    | yes              |
| Not identical to DistanceOnly               | yes              |
| Produces plausible immune biology           | yes              |
| Has interpretable behavior                  | yes              |
| Still preserves some non-promoter links     | yes, but limited |

Main caveat:

**The current model is quite proximal-biased.**

This is acceptable for the pilot, but future versions should test whether strong distal regulatory links can be rescued without losing enrichment.

---

## Current decision

Continue the project.

Do not keep tuning `lambda_distance` or `alpha_tf` right now. Current locked setting remains:

| Parameter         | Value |
| ----------------- | ----: |
| `lambda_distance` |   0.3 |
| `alpha_tf`        |   0.5 |

These settings give good enrichment while avoiding the stronger promoter collapse seen with higher distance weighting.

ArchR remains excluded from the prototype benchmark because it was operationally attempted but did not produce usable links under the current setup. SCENT is included via chromosome-wise sweep results.

---

## Next steps

### 1. Freeze current benchmark branch

Commit the current working state:

* reranker outputs
* SCENT chromosome sweep loader
* benchmark script
* ORA add-on script
* benchmark plots/tables
* developer notes

Suggested commit message:

```text
Add reranker benchmark against SCENT chromosome sweep
```

### 2. Make the repo reproducible

Turn the current benchmark branch into a clean workflow before adding more science.

Target structure:

```text
config/
workflow/
scripts/
containers/
docs/
results/
```

Use the same practical style as the previous `scRNAseq-pbmc-workflow` repo.

Minimum workflow targets:

```text
run_reranker
run_scent_chr_sweep
run_benchmark
run_ora
```

Start from 10x H5 + fragments. Do not add raw FASTQ processing yet unless it is easy.

### 3. Test transfer on breast cancer data

Run the locked reranker on the breast cancer dataset without retuning.

Use:

```text
lambda_distance = 0.3
alpha_tf = 0.5
```

Compare:

* LinkPeaks
* Coactivity
* DistanceOnly
* Reranker
* SCENT if feasible

Main question:

**Does the reranker still improve biological coherence outside the immune/PBMC context?**

### 4. Start standalone candidate generation after transfer test

Only after the repo is reproducible and breast cancer has been tested, start the standalone model.

Add:

```text
candidate_mode = "linkpeaks"
candidate_mode = "window"
```

First standalone version:

* same chromosome
* fixed cis window, initially 100–500 kb
* gene expression prevalence filter
* peak accessibility prevalence filter
* current scoring formula unchanged

Main question:

**Does the current score still work when LinkPeaks no longer defines the candidate universe?**

### 5. Later: cluster-aware / metacell version

After standalone cis-window scoring works:

* cluster cells
* optionally annotate clusters
* create metacells within clusters
* compute coactivity per cluster/metacell
* score links per cluster
* aggregate with max or best-cluster score

Do not start this until the simpler standalone version works.

---

## Current project status

The project is still a prototype, but the benchmark results are coherent enough to justify continuing.

The current model is not ready to claim validated enhancer–gene links, but it is strong enough to support the next phase:

1. reproducible repo/workflow
2. second dataset transfer test
3. standalone cis-window candidate model
