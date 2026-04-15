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


