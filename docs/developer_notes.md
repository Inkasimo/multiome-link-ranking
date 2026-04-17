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

#  AAAAAAAA

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


