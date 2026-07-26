# TODO — repository cleanup and release preparation

Scope: clean, document, and package the **current PBMC LinkPeaks-candidate peak–gene reranking
benchmark**. This is not a plan to implement the future standalone method.

Paths below are relative to the repository root. Two sources are used:

- **archive** = `multiome_test_clean_reranker_handoff_20260725_123328.tar.gz` (401 files, 94 MB unpacked)
- **tree** = `repo_tree_for_claude.txt` (1,342 entries, the live working repo)

Where a path exists in the tree but not the archive, this is stated explicitly.

---

## 0. Blockers — fix before anything else

These three items mean the repository, as it stands, cannot be run or rebuilt by a third party.
Everything else in this document is cosmetic by comparison.

### 0.1 The Snakemake feature rule calls a script that does not exist

`workflow/Snakefile` rule `linkpeaks_reranker_features` invokes:

```
Rscript scripts/run_linkpeaks_reranker.R
```

No such file exists. Present instead:

- `scripts/run_linkpeaks_reranker.R.bak` (664 lines)
- `scripts/run_linkpeaks_reranker_without_scent.R` (711 lines)

**Resolution — confirmed by CLI contract check.** Rename
`scripts/run_linkpeaks_reranker_without_scent.R` → `scripts/run_linkpeaks_reranker.R`.

Evidence for that choice:

- All 19 flags the Snakefile passes are accepted by the `_without_scent.R` option parser.
  Set difference in the direction "Snakefile sends but script rejects" is empty.
- The script constructs its outputs at lines 123–125 as
  `{run_name}_link_features.csv`, `{run_name}_baseline_links_full.csv`,
  `{run_name}_baseline_links_with_distance.csv`, and touches `--done-file`.
  The Snakefile expects exactly those four at lines 133–136.
- Every flag in `run_linkpeaks_reranker.R.bak` is also present in `_without_scent.R`,
  which additionally carries `--candidate-filter` (default `positive_score`).
  That default is what produced the recorded 15,806-candidate figure.
- `grep -ci scent` returns **0 for both scripts**. The `_without_scent` suffix records a
  historical concern, not a capability difference — SCENT was always in
  `scripts/run_scent_chr_sweep.R`. The suffix is misleading and should not be preserved.

Status: **safe**, given the above.

### 0.2 The Docker image cannot be built from the archive, and does not contain the JASPAR fix

`containers/Dockerfile` contains:

```
COPY renv.lock       /work/renv.lock
COPY renv/vendor/    /work/renv/vendor/
COPY renv/activate.R /work/renv/activate.R
COPY .Rprofile       /work/.Rprofile
```

- `renv/`, `.Rprofile` — present in tree, **absent from archive**. The live repo is fine;
  the handoff tarball is not buildable. Fix the tarball script, not the Dockerfile.
- `resources/jaspar/JASPAR2022.sqlite` — present in tree, **never `COPY`-ed into the image**.

The second point is a genuine reproducibility hole. `scripts/run_linkpeaks_reranker.R`
lines 8–48 define `seed_jaspar2022_cache()`, which seeds `BiocFileCache` with the local
`resources/jaspar/JASPAR2022.sqlite` so that the `JASPAR2022` package does not attempt a
network download from `https://jaspar2022.genereg.net/download/database/JASPAR2022.sqlite`.
This is the fix that unblocked feature generation. It is currently:

- undocumented in `README.md`
- undocumented in the Dockerfile
- dependent on a file that must be present in the bind-mounted working directory at runtime

Actions:

1. Document the requirement in `README.md` and `docs/input_output_reference.md`. **Safe.**
2. Add `resources/jaspar/JASPAR2022.sqlite.sha256` verification to the run path, or a
   `make_resources` helper that fetches and verifies. **Needs confirmation** — depends on
   whether the sqlite is redistributable (see §6).
3. Decide whether the sqlite ships in the image, in Git LFS, or as a download step.
   **Needs confirmation.**

### 0.3 The proximal-control outputs are now reproducible through the workflow — RESOLVED

Previously these outputs were produced by hand with no Snakemake rule. Now wired in:

| Component | Location |
|---|---|
| Config | `config/scent_validation_min_distance.yaml` — `min_distances: "10000,25000,50000"`, `top_n_values: "50,100,200,500"`, `high_fraction: 0.10` |
| Rule | `rule scent_validation_min_distance`, `workflow/Snakefile` |
| Wrapper | `python3 run_analysis.py run_scent_validation_min_distance` |
| Inputs | `results/pbmc/scent_validation/.done`, `scent_validation_all_ranked_methods_combined.csv` |
| Outputs | `results/pbmc/scent_validation_min_distance/{.done, scent_min_distance_method_counts.csv, scent_min_distance_topN_support_summary.csv, scent_min_distance_delta_vs_linkpeaks.csv, scent_min_distance_distance_matched_enrichment.csv}` |

Light post-processing only. **Does not re-run SCENT.** Dry-run and real run both verified.

Two small items remain, neither blocking:

1. The `.done` target is **not** in `rule all_with_scent`, so the controls are reachable only
   by explicit request. Adding it is a one-line workflow-completeness change, not a new
   analysis. **Safe.**
2. No plot is produced. See `docs/results_report.md` §6 for the recommended figure and path.
   **Safe.**
---

## 1. Keep

### Source and configuration

```
workflow/Snakefile
scripts/run_linkpeaks_reranker.R              # after rename, see §0.1
scripts/evaluate_rankings.R
scripts/run_scent_chr_sweep.R
scripts/benchmark_scent_validation.R
scripts/summarize_scent_validation_min_distance.R
run_analysis.py
config/default.yaml
config/ablations.yaml
config/scent_run.yaml
config/scent_validation.yaml
containers/Dockerfile
renv.lock
LICENSE
```

### Documentation

```
docs/benchmark_summary.md      # curated; primary source for the results report
docs/folder_structure.md       # fold into docs/input_output_reference.md, then delete
docs/developer_notes.md        # keep, but rename — see §3
```

### Results — keep in the repository

```
results/pbmc/features/pbmc_link_features.csv                    # 2.3 MB, 5,000 rows
results/pbmc/features/pbmc_baseline_links_full.csv              # 1.9 MB, 15,806 rows
results/pbmc/features/pbmc_baseline_links_with_distance.csv     # 3.5 MB, 15,806 rows
results/pbmc/scent_validation/*.csv  except the combined table  # see §4
results/pbmc/scent_validation/*.png                             # 3 plots, 330 KB
results/pbmc/scent_validation_min_distance/*.csv                # 36 KB, all four
results/pbmc/scent_chr_sweep_100kb_frac020_1000cells/scent_chr_sweep_summary.csv
```

### Results — keep per ranking mode, drop the rest (see §4)

For each of the 11 directories under `results/pbmc/rankings/`:

```
pbmc_<MODE>_summary_metrics.csv
pbmc_<MODE>_tier_summary.csv
pbmc_<MODE>_top100_links.csv
pbmc_<MODE>_distance_distribution.png
pbmc_<MODE>_validation_topN_distance_diversity_summary.csv
pbmc_<MODE>_validation_component_correlations.csv
pbmc_<MODE>_validation_distance_matched_feature_contrast.csv
```

### Scripts to keep but relocate

```
make_handoff_tarball.sh      # tree only, absent from archive → move to scripts/
make_validation_tarball.sh   # tree only, absent from archive → move to scripts/
wrapper-requirements.txt     # tree only, absent from archive → keep at root
```

---

## 2. Remove

### 2.1 Safe deletions — backup and superseded files

```
config/ablations.yaml.bak                          # 6 modes; current file has 11
workflow/Snakefile.bak                             # 3,948 B vs current 12,467 B
workflow/Snakefile_without_scent                   # superseded; SCENT is now a rule
scripts/evaluate_rankings.R.bak
scripts/evaluate_rankings_most_recent.R.bak
scripts/evaluate_rankings_recent.R.bak
scripts/run_linkpeaks_reranker.R.bak               # after §0.1 rename
run_analysis.py.bak                                # tree only
run_analysis_before_scent.py                       # tree only
renv.lock.before_scent                             # tree only
```

Status: **safe.** All are strict subsets or predecessors of retained files, verified by
flag-set comparison for the R scripts.

### 2.2 Safe deletions — empty files

```
config/config.yaml            # 0 bytes. Note: run_analysis.py defaults to
                              # --configfile config/default.yaml, so this is dead.
containers/environment.yml    # 0 bytes. The stack is renv + Dockerfile; conda is unused.
```

Status: **safe.** Confirm no CI or documentation references `config/config.yaml` first.

### 2.3 Safe deletions — scratch and inventory artifacts

```
scent_sweep_file_inventory.txt      # tree only
repo_tree_for_claude.txt            # tree only; a handoff aid, not a repo artifact
scripts/run_analysis.py             # byte-identical to root run_analysis.py (diff confirms)
```

Status: **safe.** Keep the root copy; `run_analysis.py` is documented as a root-level entry point.

### 2.4 Large deletions — nested and backup result trees

```
tarballs/                                            # 630 tree entries
  tarballs/multiome_test_clean_reranker_handoff_20260725_123328.tar.gz
  tarballs/pbmc_ablation_plus_lambda_results_check.tar.gz
  tarballs/pbmc_ablation_results_check.tar.gz
  tarballs/pbmc_reranker_validation_20260629_185949.tar.gz
  tarballs/pbmc_reranker_validation_20260629_190106.tar.gz
  tarballs/pbmc_reranker_validation_outputs.tar.gz
  tarballs/pbmc_scent_validation_review.tar.gz
  tarballs/results/                                  # full nested copy of results/pbmc
results/pbmc/rankings_backup/                        # 7-mode copy; predates moddist modes
results/pbmc.before_restore_20260715_195026/         # empty except one subdirectory
```

`tarballs/` alone accounts for roughly half the entries in the repository tree, and
`tarballs/results/` duplicates `results/pbmc/` wholesale. `rankings_backup/` contains only
`coactivity, coactivity_distance, coactivity_tf, distance_only, full, full_lambda_0_1,
full_lambda_0_2, linkpeaks` — it lacks the `full_moddist_*` and `distance_mod_only_*`
directories, confirming it predates the modified-distance work and is not a fallback.

Status: **needs confirmation** for `tarballs/*.tar.gz` — these are the only copies of some
intermediate states. Archive them off-repo (see §5) before deleting.
Status: **safe** for `tarballs/results/`, `rankings_backup/`, `pbmc.before_restore_*`.

Add to `.gitignore`:

```
tarballs/
results/**/rankings_backup/
results/**/*.before_restore_*/
```

### 2.5 Documentation to remove or convert

```
docs/What this pipeline does.docx     # 24 KB, filename contains spaces
docs/math.docx                        # 25 KB
docs/plan.docx                        # 30 KB
```

Binary Office files in a documentation directory cannot be diffed, reviewed, or read on
GitHub. The filename `What this pipeline does.docx` will break naive shell tooling.

Action: extract to text, merge any unique content into `docs/method_report.md` or
`docs/future_standalone_v0.md`, then delete. **Needs confirmation** — I could not read
these during preparation, so I cannot confirm they contain nothing unique. Check before
deleting, especially `math.docx` against `docs/method_report.md`.

```
README_SPLIT.md
```

Stale. It documents `scripts/run_linkpeaks_reranker.R` and `scripts/evaluate_rankings.R`
as a "proposed split" — that split is now the implemented architecture, so the document
describes the present as a proposal. Its Docker and CLI examples are the only content worth
keeping and they belong in `README.md`. Status: **safe** after §5 step 5.

---

## 3. Rename and move

| From | To | Rationale | Status |
|---|---|---|---|
| `scripts/run_linkpeaks_reranker_without_scent.R` | `scripts/run_linkpeaks_reranker.R` | Matches Snakefile; suffix is a misnomer (§0.1) | safe |
| `docs/developer_notes.md` | `docs/lab_notebook.md` | 66 KB of unedited working notes containing "not publishable as-is", "Stop polishing this reranker", "Option A — Archive / stop". Valuable and honest, but must not read as the project's conclusion. Add a header stating it is an unedited chronological notebook. | safe |
| `config/scent_validation_with_producer.yaml` | `config/scent_validation.yaml` | The `_with_producer` file is the documented superset (carries `scent_support_rule`, `scent_p_threshold`, `scent_min_score` and explanatory comments). The current `scent_validation.yaml` is a reduced, machine-written form that **omits those three keys** — they only survive via Snakefile defaults at lines 122–123. Promote the documented file. | needs confirmation |
| `make_handoff_tarball.sh` | `scripts/make_handoff_tarball.sh` | Root clutter | safe |
| `make_validation_tarball.sh` | `scripts/make_validation_tarball.sh` | Root clutter | safe |
| `docs/folder_structure.md` | merged into `docs/input_output_reference.md` | Overlapping scope; the standalone file describes the layout as "proposed" | safe |

---

## 4. Stale and duplicated outputs to delete

### 4.1 Duplicated LinkPeaks baseline artifacts — the largest single win

Six file patterns are **byte-identical across all 11 ranking directories** (md5-verified).
Each mode's evaluation re-emits the same baseline, because the baseline does not depend on
the score mode.

| Pattern | Size each | Copies | Redundant |
|---|---|---|---|
| `pbmc_<MODE>_linkpeaks_baseline_distance_same_universe.csv` | 1,115,427 B | 11 | 10.6 MB |
| `pbmc_<MODE>_linkpeaks_baseline_same_universe.csv` | 585,581 B | 11 | 5.6 MB |
| `pbmc_<MODE>_linkpeaks_baseline_ora_dotplot_same_universe.png` | 272,382 B | 11 | 2.6 MB |
| `pbmc_<MODE>_linkpeaks_gene_rank_same_universe.csv` | 75,889 B | 12 | 0.8 MB |
| `pbmc_<MODE>_linkpeaks_baseline_ora_GO_BP_same_universe.csv` | 5,579 B | 12 | 61 KB |
| `pbmc_<MODE>_validation_manifest.csv` | 738 B | 11 | 7 KB |

Total redundant: **~20.6 MB of the 62 MB `results/pbmc/rankings/` tree.**

The copy count is 12 rather than 11 for two patterns because in the `linkpeaks/` mode
directory the "model" output *is* the baseline, so `pbmc_linkpeaks_model_gene_rank.csv`
and `pbmc_linkpeaks_ora_GO_BP.csv` are additional identical copies.

Action: promote one copy of each to a shared location and delete the rest.

```
results/pbmc/baseline_same_universe/
  pbmc_linkpeaks_baseline_same_universe.csv
  pbmc_linkpeaks_baseline_distance_same_universe.csv
  pbmc_linkpeaks_gene_rank_same_universe.csv
  pbmc_linkpeaks_baseline_ora_GO_BP_same_universe.csv
  pbmc_linkpeaks_baseline_ora_dotplot_same_universe.png
results/pbmc/rankings/validation_manifest.csv
```

Status: **safe** for deletion of the redundant copies.
Status: **needs confirmation** for changing `scripts/evaluate_rankings.R` to stop
re-emitting them — that is a code change and will alter future output layout. Recommend
doing the deletion now and filing the code change as a follow-up issue, so the committed
results and the code that generated them stay consistent for this release.

### 4.2 Per-mode outputs to drop from the repository

Per ranking directory, these are large, derivable from retained files, or of low review value:

```
pbmc_<MODE>_ranked_links.csv                              # 2.5–2.7 MB × 11 = ~29 MB
pbmc_<MODE>_validation_promoted_demoted_inspection.csv    # ~61 KB × 11
pbmc_<MODE>_top_promoted_vs_linkpeaks.csv                 # ~52 KB × 11
pbmc_<MODE>_top_demoted_vs_linkpeaks.csv                  # ~53 KB × 11
pbmc_<MODE>_model_gene_rank.csv                           # ~76 KB × 11
pbmc_<MODE>_ora_dotplot.png                               # 80–272 KB × 11
pbmc_<MODE>_linkpeaks_vs_model_scatter.png                # 131–776 KB × 11
pbmc_<MODE>_validation_distance_bin_rank_summary.csv
pbmc_<MODE>_validation_tf_motif_support_summary.csv
pbmc_<MODE>_ora_GO_BP.csv
pbmc_<MODE>_topN_overlap_vs_linkpeaks_same_universe.csv   # 42–48 bytes; a single number
```

The `*_ranked_links.csv` files are the largest item at ~29 MB total. Every one is fully
regenerable from `results/pbmc/features/pbmc_link_features.csv` plus
`config/ablations.yaml`, in minutes, without re-running the heavy Seurat step. They belong
in the external archive, not in Git.

Note on `linkpeaks/`: three files in that directory are byte-identical to each other
(`pbmc_linkpeaks_top100_links.csv`, `pbmc_linkpeaks_top_promoted_vs_linkpeaks.csv`,
`pbmc_linkpeaks_top_demoted_vs_linkpeaks.csv`, all 51,420 B) because the baseline cannot
be promoted or demoted relative to itself. Keep only `top100_links`.

Status: **safe**, provided §5 step 3 (external archive) completes first.

### 4.3 SCENT sweep outputs

```
results/pbmc/scent_chr_sweep_100kb_frac020_1000cells/chr*/scent_candidates_chr*.csv   # 22 files
results/pbmc/scent_chr_sweep_100kb_frac020_1000cells/chr*/scent_links_chr*.csv        # 22 files
results/pbmc/scent_chr_sweep_100kb_frac020_1000cells/scent_links_all_chromosomes.csv  # 6.7 MB
```

17 MB total. Keep `scent_chr_sweep_summary.csv` (per-chromosome status, row counts, score
ranges) in the repository; move the rest to the external archive.

**Do not discard the candidate files.** The 22 `scent_candidates_chr*.csv` files together
form a 117,811-pair de novo cis-window candidate universe (9,891 genes, 46,936 peaks)
generated by `scripts/run_scent_chr_sweep.R` lines 325–365 — independent of LinkPeaks. They
are the starting point for the next phase and are described in
`docs/future_standalone_v0.md`. Archive them; do not delete them.

Status: **safe** to move out of Git. **Do not delete outright.**

### 4.4 Large validation table

```
results/pbmc/scent_validation/scent_validation_all_ranked_methods_combined.csv   # 6.1 MB
```

This is the input to `summarize_scent_validation_min_distance.R`, so it is not purely
derived output — it is a pipeline intermediate with a downstream consumer. Move to the
external archive and document the dependency in `docs/input_output_reference.md`.

Status: **needs confirmation** — removing it from Git breaks the ability to re-run the
min-distance summary from a fresh clone. Either keep it in Git LFS, or make the new
min-distance rule depend on the SCENT-validation rule rather than on the committed file.

### 4.5 Sentinel files

`.done` files exist at 14 locations and are byte-empty (`d41d8cd9...`, 11 copies under
`rankings/` alone). They are legitimate Snakemake sentinels and should be **kept** where
the corresponding outputs are kept, and gitignored where outputs move to the archive.
Do not delete them from working directories — the workflow depends on them.

### 4.6 Missing from the archive, cannot be assessed

```
results/pbmc/audits/reranker_output_audit_checks.csv
```

Present in the tree, excluded from the handoff tarball. Its columns and content could not
be inspected. Review it manually and decide keep/remove; if it contains automated
consistency checks it should be kept and documented in
`docs/input_output_reference.md`. **Needs confirmation.**

---

## 5. Exact order of operations

Do these in sequence. Steps 1–4 are non-destructive.

**1. Tag the current state.**

```bash
git tag pre-cleanup-20260725
git push --tags
```

Nothing below is recoverable without this if something is deleted before being archived.

**2. Fix the blocker rename and verify the workflow resolves.**

```bash
git mv scripts/run_linkpeaks_reranker_without_scent.R scripts/run_linkpeaks_reranker.R
python3 run_analysis.py run_all_score_modes --dry-run
```

A clean dry run confirms §0.1. Do not proceed until it passes.

**3. Build the external archive.**

Copy, do not move, the following to a staging directory for Zenodo upload:

```
results/pbmc/rankings/*/pbmc_*_ranked_links.csv
results/pbmc/rankings/*/pbmc_*_ora_dotplot.png
results/pbmc/rankings/*/pbmc_*_linkpeaks_vs_model_scatter.png
results/pbmc/rankings/*/pbmc_*_validation_promoted_demoted_inspection.csv
results/pbmc/rankings/*/pbmc_*_top_{promoted,demoted}_vs_linkpeaks.csv
results/pbmc/scent_chr_sweep_100kb_frac020_1000cells/
results/pbmc/scent_validation/scent_validation_all_ranked_methods_combined.csv
tarballs/*.tar.gz
```

Record `sha256sum` for every file. Nothing is deleted until the archive is uploaded and
the DOI is minted (see `docs/release_checklist.md`).

**4. Write `.gitignore` before deleting anything.**

Otherwise the deletions get re-added on the next `git add -A`.

```
data/
resources/jaspar/*.sqlite
tarballs/
renv/library/
renv/staging/
results/**/rankings_backup/
results/**/*.before_restore_*/
results/pbmc/rankings/*/pbmc_*_ranked_links.csv
results/pbmc/scent_chr_sweep_*/chr*/
results/pbmc/scent_validation/scent_validation_all_ranked_methods_combined.csv
*.bak
```

**5. Safe deletions.** §2.1, §2.2, §2.3, and the `tarballs/results/`,
`rankings_backup/`, `pbmc.before_restore_*` items from §2.4. Then delete `README_SPLIT.md`
after its Docker and CLI examples have been merged into `README.md`.

**6. Renames and moves.** §3, excluding the two items marked *needs confirmation*.

**7. Deduplicate the baseline artifacts.** §4.1. Promote one copy of each of the six
patterns, delete the other 10–11, commit as a single reviewable change.

**8. Drop archived per-mode outputs.** §4.2, §4.3, §4.4 — only after step 3 has a DOI.

**9. Add the missing documentation.** §7 below.

**10. Add reproducibility metadata.** §8 below.

**11. Confirm-first items.** Work through everything marked *needs confirmation*, including
the three-key `scent_validation.yaml` promotion and the `.docx` extraction.

**12. Verify from a clean clone.**

```bash
git clone <repo> /tmp/verify && cd /tmp/verify
docker build -t multiome-reranker:test -f containers/Dockerfile .
python3 run_analysis.py run_all_score_modes --dry-run
```

Both must succeed with no files outside the clone, other than `data/` and the JASPAR sqlite.

---

## 6. Results that should be summarised rather than committed

These are already summarised in the documents produced alongside this file. Each raw source
should go to the external archive with a pointer left in the docs.

| Raw source | Summarised in |
|---|---|
| `results/pbmc/rankings/*/pbmc_*_ranked_links.csv` (11 × ~2.6 MB) | `docs/results_report.md` — ranking-mode inventory |
| `results/pbmc/scent_chr_sweep_*/chr*/scent_links_chr*.csv` (22 files) | `scent_chr_sweep_summary.csv` + `docs/results_report.md` |
| `results/pbmc/scent_validation/scent_validation_all_ranked_methods_combined.csv` (6.1 MB) | `scent_validation_topN_support_summary.csv`, `scent_validation_method_counts.csv` |
| `results/pbmc/rankings/*/pbmc_*_validation_promoted_demoted_inspection.csv` | `docs/results_report.md` — qualitative note only |
| `docs/benchmark_summary.md` (29 KB, chronological, contains superseded verdicts) | `docs/results_report.md` — keep the original as the record |

---

## 7. Missing documentation to write

| File | Status |
|---|---|
| `README.md` | rewrite — current version describes the pipeline as `scripts/multiome_rie_v1.R`, a file that now lives in `scripts/dev/` and is not part of the workflow |
| `docs/method_report.md` | new |
| `docs/results_report.md` | new |
| `docs/input_output_reference.md` | new |
| `docs/release_checklist.md` | new |
| `docs/future_standalone_v0.md` | new |
| `docs/competitor_positioning.md` | new |
| `CITATION.cff` | new — absent from both tree and archive |
| `docs/lab_notebook.md` header | new — one paragraph stating it is unedited notes |
| `.gitattributes` | new, if Git LFS is used |

The existing `README.md` is actively misleading in three places and must not ship as-is:

1. It documents a run command against `scripts/multiome_rie_v1.R`. That file is in
   `scripts/dev/` and is not invoked by the workflow.
2. Its "Expected outputs" list uses a flat `results/*_ranked_links.csv` layout. The actual
   layout is `results/<dataset>/rankings/<mode>/`.
3. It presents the scoring formula with `λ = 0.30, α = 0.50` as the default. `config/default.yaml`
   does set those, but every result carried forward in the SCENT validation uses `λ = 0.10`.

---

## 8. Missing reproducibility metadata

| Item | Determinable from archive? | Action |
|---|---|---|
| Genome build | **Yes** — hg38, via `EnsDb.Hsapiens.v86` and `BSgenome.Hsapiens.UCSC.hg38` (`run_linkpeaks_reranker.R` lines 375–384, 446) | state in README and method report |
| Motif source | **Yes** — JASPAR2022, `collection=CORE`, `species=9606`, `tax_group=vertebrates` (line 616) | state in method report |
| R and Bioconductor | **Yes** — R 4.4.3, Bioconductor 3.20 (`renv.lock`, Dockerfile) | version table in release checklist |
| Key package versions | **Yes** — 268 packages pinned in `renv.lock` | version table in release checklist |
| SCENT provenance | **Yes** — `immunogenomics/SCENT`, v1.0.1, commit `e80b5ba6b445f972c7fe28fb41e24ef4f5b2e373` | cite exactly, not "SCENT (GitHub)" |
| Random seed | **Yes** — `seed: 42` in `config/default.yaml`, threaded to all four scripts | state in README |
| **Dataset accession** | **Yes, recovered** — `pbmc_unsorted_10k`, Cell Ranger ARC 2.0.0, reference `GRCh38-2020-A`, from the `atac_fragments.tsv.gz` header | record in `config/default.yaml`; see `dataset_provenance.yaml` |
| **Cell count** | **Yes, by inference** — there is no QC step, so N = barcodes in the filtered matrix ≈ 12,012 | confirm with `ncol()` on a rerun |
| Cell-level QC | **None performed** — no `min.cells`, `min.features` or `subset()` | document as a limitation |
| Annotation consistency | **Mismatched** — counts from `GRCh38-2020-A`, TSS from `EnsDb.Hsapiens.v86` | document as a limitation |
| `--candidate-filter` value used | Inferable — default `positive_score`, consistent with 15,806 | promote to `config/default.yaml` so it is explicit |
| `--min-distances`, `--top-n-values`, `--high-fraction` | **Yes, resolved** — version-controlled in `config/scent_validation_min_distance.yaml` | no action |
| Docker image digest | **No** — no image was published | build, push, record digest |

---

## 9. Cannot be determined from the archive

These block or constrain the release and require your input.

### Resolved since first draft

1. ~~**Dataset accession or source.**~~ **RESOLVED.** Recovered from the `cellranger-arc` header
   of `data/atac_fragments.tsv.gz`:

   | | |
   |---|---|
   | Sample ID | `pbmc_unsorted_10k` |
   | Title | PBMC from a Healthy Donor - No Cell Sorting (10k) |
   | Pipeline | `cellranger-arc count`, Cell Ranger ARC 2.0.0 |
   | Reference | `refdata-cellranger-arc-GRCh38-2020-A-2.0.0`, version `2020-A` |
   | `reference_fasta_hash` | `b6f131840f9f337e7b858c3d1e89d7ce0321b243` |
   | `reference_gtf_hash` | `3b4c36ca3bade222a5b53394e8c07a18db7ebb11` |
   | Dataset page | https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-no-cell-sorting-10-k-1-standard-2-0-0 |
   | Download base | `https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_unsorted_10k/` |
   | Estimated cells | 12,012 |
   | ATAC peaks called | 111,857 (50,918 recoverable from results = 45.5%) |
   | Published | 2021-05-03 |
   | **License** | **CC BY 4.0** |

   Note this is **not** `pbmc_granulocyte_sorted_10k` from the Signac multiome vignette, despite
   identical local filenames. Fragment-file sizes differ (2,917,757,251 B versus 2,051,027,831 B).

   - [ ] Add `dataset_provenance.yaml` content to `config/default.yaml` as comments
   - [ ] Run `sha256sum data/*` **now**, while the files are on disk

3. ~~**Whether `data/` is redistributable.**~~ **RESOLVED — CC BY 4.0.** Derived outputs may be
   redistributed with attribution. Recommendation is still not to deposit the raw data, since 10x
   hosts it durably; cite the accession and pin checksums instead.

**Still blocking:**

2. **JASPAR2022 sqlite redistribution terms.** Determines whether the file can ship in the
   repo, in LFS, in the Docker image, or must be fetched at runtime.

**Needed for release metadata:**

4. Author name as it should appear, ORCID, affiliation — for `CITATION.cff`. The only
   attribution in the archive is `Copyright (c) 2026 Inkasimo` in `LICENSE` and
   `LABEL maintainer="Inkasimo"` in the Dockerfile.
5. Final repository name. Recommended: `multiome-peak-gene-reranking-benchmark` — names it
   as a benchmark rather than a method.
6. Final Docker image name and registry. By analogy with the existing
   `ghcr.io/inkasimo/scrnaseq-pbmc-workflow:v2.0.0`, suggest
   `ghcr.io/inkasimo/multiome-peak-gene-reranking-benchmark:v0.1.0`.
   The Dockerfile currently carries `org.opencontainers.image.title="multiome-reranking-benchmark"`
   and `run_analysis.py` defaults to `DEFAULT_IMAGE = "multiome-reranking-benchmark:v0.1.0"` — these
   must be reconciled with whatever is chosen.
7. Git release tag. Suggest `v0.1.0`, reflecting benchmark-not-method status.
8. Zenodo concept DOI target.
9. **LICENSE for results and documentation.** The code is MIT. Data and figure reuse terms
   are unstated; CC-BY-4.0 for `docs/` and `results/` is the conventional pairing.

**Governs archival decisions:**

10. Whether the ~50 MB of retained results is acceptable in Git, or whether the threshold
    should be lower. After §4 the repository is roughly 12–15 MB of results; before §4 it
    is 94 MB plus the `tarballs/` tree.
11. Whether Git LFS is available and wanted, specifically for
    `scent_validation_all_ranked_methods_combined.csv` (6.1 MB) and
    `resources/jaspar/JASPAR2022.sqlite`.
12. Whether `docs/lab_notebook.md` (formerly `developer_notes.md`) should be public at all.
    It contains candid negative self-assessment. My recommendation is to publish it with a
    header — it is the most honest artifact in the repository — but it is your call.

**External citation targets to fix:**

13. Exact citations for Signac/LinkPeaks, SCENT, ArchR, Cicero, SCARlink, CREMA, TRIPOD,
    SCENIC+, Pando, LINGER, FigR, JASPAR2022, EnsDb, clusterProfiler. `docs/competitor_positioning.md`
    names each method but deliberately does not invent bibliographic details.

**New findings from provenance recovery — document, do not fix in this release:**

16. **No cell-level QC is performed.** `run_linkpeaks_reranker.R:380` calls
    `CreateSeuratObject()` without `min.cells` or `min.features`, and no `subset()` step exists
    anywhere. All ~12,012 barcodes enter the analysis. Comparable Signac multiome workflows filter
    on `nCount_ATAC`, `nCount_RNA`, `nucleosome_signal` and `TSS.enrichment` first. Documented in
    `docs/method_report.md` §2 and §14. **Do not add QC now** — it would invalidate every
    committed result. File as a v0.2 issue.
17. **Annotation versions are mismatched.** Counts come from the Cell Ranger ARC `GRCh38-2020-A`
    reference; TSS coordinates driving `distance_bp` come from `EnsDb.Hsapiens.v86` (Ensembl 86,
    2016). Genes absent from Ensembl 86 are silently dropped and TSS positions may drift. Since
    distance is the central confound in this benchmark, state it explicitly. Confirm the exact
    GENCODE/Ensembl correspondence of `2020-A` from 10x's reference build notes before citing it.

**Content checks:**

14. Whether `docs/What this pipeline does.docx`, `docs/math.docx`, `docs/plan.docx` contain
    anything not already in the Markdown docs.
15. Contents and purpose of `results/pbmc/audits/reranker_output_audit_checks.csv`.
