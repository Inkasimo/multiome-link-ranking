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

