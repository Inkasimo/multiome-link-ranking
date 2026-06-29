
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
