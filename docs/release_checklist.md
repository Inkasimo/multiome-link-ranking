# Release checklist

Release-readiness for the PBMC LinkPeaks-candidate peak–gene reranking benchmark.

Complete `TODO.md` §0 (the three blockers) and §5 (ordered cleanup) before starting here.
Items marked **BLOCKING** must be resolved before the repository is made public.

---

## 1. License

### Code — already in place

`LICENSE` is MIT, `Copyright (c) 2026 Inkasimo`. Consistent with the existing
`Inkasimo/scRNAseq-pbmc-workflow` repository. No change needed.

- [ ] Confirm the copyright holder name is how it should appear publicly

### Documentation and results — not yet chosen

MIT is a software license and is a poor fit for prose and data tables. Recommended:

- Code (`scripts/`, `workflow/`, `containers/`, `run_analysis.py`, `config/`) → **MIT**, as now
- Documentation (`docs/`, `README.md`) and results (`results/`) → **CC-BY-4.0**

- [ ] Add a `## License` section to `README.md` stating the split
- [ ] Add `LICENSE-docs` containing CC-BY-4.0, or state the split in `README.md` alone

### Third-party — BLOCKING

- [ ] **`resources/jaspar/JASPAR2022.sqlite` redistribution terms.** Determines whether the file
      can ship in the repository, in Git LFS, inside the Docker image, or must be fetched at
      runtime. Until resolved, keep it gitignored.
- [ ] **`data/` redistribution terms.** Governs the Zenodo deposit contents.
- [ ] Confirm SCENT's license permits redistributing derived outputs
      (`immunogenomics/SCENT`, commit `e80b5ba6b445f972c7fe28fb41e24ef4f5b2e373`)

---

## 2. CITATION.cff

Not present in either the archive or the repository tree. Create at the root.

```yaml
cff-version: 1.2.0
message: "If you use this benchmark, please cite it as below."
type: software
title: "<FINAL_REPO_NAME>"
abstract: >
  A reproducible, containerized benchmark that reranks Signac LinkPeaks candidate
  peak-gene pairs from 10x PBMC multiome data using interpretable RNA-ATAC coactivity,
  genomic distance and TF/motif support, and evaluates the result against SCENT with
  explicit controls for promoter-proximity confounding.
authors:
  - family-names: "<FAMILY_NAME>"        # TO FILL
    given-names: "<GIVEN_NAME>"          # TO FILL
    orcid: "https://orcid.org/<ORCID>"   # TO FILL
    affiliation: "<AFFILIATION>"         # TO FILL
version: "0.1.0"
date-released: "<YYYY-MM-DD>"            # TO FILL
license: MIT
repository-code: "https://github.com/<USER>/<REPO>"
doi: "<ZENODO_DOI>"                      # TO FILL after step 3
keywords:
  - single-cell multiome
  - peak-gene links
  - scATAC-seq
  - scRNA-seq
  - benchmark
  - reproducible research
  - Snakemake
```

- [ ] Author family name, given name — **the only attribution anywhere in the archive is
      `Inkasimo`**, in `LICENSE` and the Dockerfile `LABEL maintainer`
- [ ] ORCID
- [ ] Affiliation
- [ ] `date-released`
- [ ] Final repository name and URL
- [ ] Zenodo DOI, after step 3
- [ ] Validate: `cffconvert --validate`

The `abstract` deliberately says "benchmark" and "evaluates", not "method" or "improves". Keep it
that way — `CITATION.cff` is machine-readable and gets scraped.

---

## 3. Zenodo DOI

- [ ] Decide the deposit scope. Recommended: pipeline outputs excluded from Git, **not** raw input
      data (pending §1 licensing)
- [ ] Assemble from `TODO.md` §5 step 3:
      - all 11 `results/pbmc/rankings/*/pbmc_*_ranked_links.csv` (~29 MB)
      - `results/pbmc/scent_chr_sweep_100kb_frac020_1000cells/` complete (~17 MB)
      - `results/pbmc/scent_validation/scent_validation_all_ranked_methods_combined.csv` (6.1 MB)
      - per-mode `ora_dotplot.png` and `linkpeaks_vs_model_scatter.png` (~7 MB)
      - per-mode `validation_promoted_demoted_inspection.csv`, `top_{promoted,demoted}` tables
      - `tarballs/*.tar.gz`, if they are to be preserved at all
- [ ] **Include the 22 `scent_candidates_chr*.csv` files.** These constitute the 117,811-pair
      de novo cis-window candidate universe and are the starting point for the next phase
      (`docs/future_standalone_v0.md` §3). Do not lose them
- [ ] Generate `sha256sum` for every file; ship a `MANIFEST.sha256`
- [ ] Include a `README` in the deposit describing the tree and pointing at the GitHub repository
- [ ] Enable the GitHub–Zenodo integration on the repository **before** creating the release tag,
      so the tag triggers archival automatically
- [ ] Reserve the DOI, insert into `CITATION.cff` and the `README.md` badge, commit, then publish
- [ ] Record the **concept DOI** (version-independent) in `README.md`, not only the version DOI

---

## 4. Docker image name and versioning

Three names currently disagree and must be reconciled:

| Location | Current value |
|---|---|
| `containers/Dockerfile` | `org.opencontainers.image.title="multiome-link-ranking"` |
| `run_analysis.py` | `DEFAULT_IMAGE = "multiome-link-ranking:pilot"` |
| `README_SPLIT.md` | `multiome-link-ranking:pilot` |

Recommended, by analogy with `ghcr.io/inkasimo/scrnaseq-pbmc-workflow:v2.0.0`:

```
ghcr.io/inkasimo/multiome-peak-gene-reranking-benchmark:v0.1.0
```

- [ ] Choose the final name — **BLOCKING** for §5 and §6
- [ ] Update `containers/Dockerfile` labels: `title`, `description`, `version`, `source`, `licenses`
- [ ] Update `DEFAULT_IMAGE` in `run_analysis.py`
- [ ] Update all `README.md` examples
- [ ] Tag `v0.1.0` and `latest`
- [ ] Push to GHCR; make the package public
- [ ] **Record the image digest** (`sha256:...`) in `README.md`. A tag is mutable; a digest is not
- [ ] Verify the image builds from a clean clone — note that the Dockerfile `COPY`s `renv/vendor/`,
      `renv/activate.R` and `.Rprofile`, which were absent from the handoff tarball (`TODO.md` §0.2)
- [ ] Decide whether `resources/jaspar/JASPAR2022.sqlite` ships in the image. It currently does
      not, and feature generation requires it in the bind-mounted working directory

Version `0.1.0`, not `1.0.0`: this is a frozen benchmark, not a released method.

---

## 5. GitHub release

- [ ] Repository name — recommended `multiome-peak-gene-reranking-benchmark`. The name is the
      claim; it should say "benchmark" without anyone having to read the README
- [ ] Repository description, one line, e.g. "Reproducible benchmark reranking LinkPeaks
      candidate peak–gene links from PBMC multiome data, with SCENT validation and explicit
      promoter-proximity controls."
- [ ] Topics: `single-cell`, `multiome`, `scatac-seq`, `scrna-seq`, `peak-gene-links`,
      `benchmark`, `snakemake`, `reproducible-research`, `bioinformatics`, `seurat`, `signac`
- [ ] Tag `v0.1.0` after cleanup is committed
- [ ] Release title: "v0.1.0 — frozen reranking benchmark"
- [ ] Release notes must state: benchmark not method; headline result with the distance-only
      caveat; the 100 kb validator-window limitation; link to `docs/results_report.md` §11
      ("What cannot be claimed")
- [ ] Attach nothing large to the release; point at Zenodo
- [ ] Confirm Zenodo archived the tag and the DOI resolves

---

## 6. README badges

Add once the corresponding step is complete. Do not add a badge that does not resolve.

```markdown
[![DOI](https://zenodo.org/badge/DOI/<DOI>.svg)](https://doi.org/<DOI>)
[![Snakemake](https://img.shields.io/badge/snakemake-7.32.4-brightgreen.svg)](https://snakemake.github.io)
[![Docker](https://img.shields.io/badge/docker-ghcr.io-blue.svg)](<PACKAGE_URL>)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-4.4.3-blue.svg)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.20-blue.svg)](https://bioconductor.org/)
```

- [ ] DOI — after §3
- [ ] Snakemake 7.32.4 — verified from the Dockerfile pip install
- [ ] Docker/GHCR — after §4
- [ ] License MIT
- [ ] R 4.4.3, Bioconductor 3.20 — verified from `renv.lock`
- [ ] **No CI badge unless CI is actually added.** There is no `.github/workflows/` in this
      repository. Consider a dry-run-only workflow (`snakemake --dry-run` against a stub config),
      which is cheap and catches the class of breakage described in `TODO.md` §0.1

---

## 7. Environment and package versions

All values below verified from `renv.lock` (268 packages) and `containers/Dockerfile`.

### Runtime

| Component | Version | Source |
|---|---|---|
| Base image | `rocker/r-ver:4.4.3` | Dockerfile |
| R | 4.4.3 | `renv.lock` |
| Bioconductor | 3.20 | `renv.lock`, pinned via `renv::settings$bioconductor.version()` |
| Snakemake | 7.32.4 | Dockerfile pip |
| PuLP | 2.7.0 | Dockerfile pip |
| PyYAML | ≥ 6.0 | Dockerfile pip |
| Python | 3 (system, Ubuntu), venv at `/opt/venv` | Dockerfile |

### Key R packages

| Package | Version | Source |
|---|---|---|
| Seurat | 5.4.0 | vendored (`renv/vendor/Seurat_5.4.0.tar.gz`) |
| SeuratObject | 5.4.0 | vendored |
| Signac | 1.16.0 | vendored (`renv/vendor/Signac_1.16.0.tar.gz`) |
| TFMPvalue | 0.0.9 | vendored |
| SCENT | 1.0.1 | GitHub `immunogenomics/SCENT`, commit `e80b5ba6b445f972c7fe28fb41e24ef4f5b2e373` |
| JASPAR2022 | 0.99.8 | Bioconductor |
| TFBSTools | 1.44.0 | Bioconductor |
| motifmatchr | 1.28.0 | Bioconductor |
| clusterProfiler | 4.14.6 | Bioconductor |
| org.Hs.eg.db | 3.20.0 | Bioconductor |
| EnsDb.Hsapiens.v86 | 2.99.0 | Bioconductor |
| BSgenome.Hsapiens.UCSC.hg38 | 1.4.5 | Bioconductor |
| GenomicRanges | 1.58.0 | Bioconductor |
| Matrix | 1.7-5 | CRAN |
| data.table | 1.18.4 | CRAN |
| lme4 | 2.0-1 | CRAN |
| RSQLite | 3.53.2 | CRAN |
| optparse | 1.8.2 | CRAN |
| ggplot2 | 3.5.2 | CRAN |

### System tools

`tabix`, `samtools`, `bedtools`, plus `libhdf5-dev`, `libglpk-dev`, `libgsl-dev`, `libnlopt-dev`
and the usual graphics stack. Full list in `containers/Dockerfile`.

- [ ] Add the tables above to `README.md` or a `docs/environment.md`
- [ ] Note that four packages are **vendored** as local tarballs under `renv/vendor/` and that
      the Dockerfile `COPY`s that directory — a clone without it cannot build
- [ ] Record `renv.lock`'s own `sha256sum` in the release notes

---

## 8. Reproducibility checklist

- [ ] `git tag pre-cleanup-20260725` exists and is pushed
- [ ] **BLOCKING — dataset accession or source recorded.** Nothing in `config/`, `scripts/`,
      `docs/benchmark_summary.md` or `docs/developer_notes.md` identifies which PBMC multiome
      dataset `filtered_feature_bc_matrix.h5` and `atac_fragments.tsv.gz` came from. A public
      release without this is not reproducible
- [ ] Genome build stated in `README.md` — hg38, `EnsDb.Hsapiens.v86`,
      `BSgenome.Hsapiens.UCSC.hg38`
- [ ] Motif source stated — JASPAR2022 `CORE`, `tax_group=vertebrates`, `species=9606`
- [ ] Seed stated — `seed: 42` in `config/default.yaml`, threaded to all four scripts
- [ ] JASPAR sqlite requirement documented, including that it must be in the **bind-mounted**
      working directory (`TODO.md` §0.2)
- [ ] `resources/jaspar/JASPAR2022.sqlite.sha256` present and verified
- [ ] `--candidate-filter` promoted to `config/default.yaml`. It defines the candidate universe
      (15,806 pairs at `positive_score`) and is currently an invisible script default
- [ ] `run_analysis.py <section> --dry-run` succeeds for every section from a clean clone
- [ ] `docker build` succeeds from a clean clone
- [x] **Snakemake rule added for the min-distance controls — DONE.**
      `rule scent_validation_min_distance` (`workflow/Snakefile`), configured by
      `config/scent_validation_min_distance.yaml`, exposed as
      `python3 run_analysis.py run_scent_validation_min_distance`. Dry-run and real run passed
- [ ] Optional: add the min-distance `.done` target to `rule all_with_scent`. It is currently
      reachable only by explicit request, so `all_with_scent` does not regenerate it
- [x] **`--min-distances` / `--high-fraction` — RESOLVED.** Now version-controlled in
      `config/scent_validation_min_distance.yaml`: `min_distances: "10000,25000,50000"`,
      `top_n_values: "50,100,200,500"`, `high_fraction: 0.10`. Outputs regenerated through the
      workflow, so the committed values and the config agree
- [ ] Cell count after QC recorded — currently nowhere in the archive
- [ ] Total runtime recorded. `scent_chr_sweep_summary.csv` has `status = skipped_existing` and
      blank `runtime_minutes` for all 22 chromosomes, so compute cost is unknown
- [ ] `wrapper-requirements.txt` present at the root (tree only, absent from the archive)
- [ ] `.gitignore` written **before** any deletion (`TODO.md` §5 step 4)
- [ ] `.gitattributes` present if Git LFS is used
- [ ] Known interface hazards documented — peak ID format mismatch, 4,976 vs 5,000, candidate
      window vs validator window, no cell-type stratification
      (`docs/input_output_reference.md` §9)

---

## 9. Result archival checklist

- [ ] Every file in the Zenodo deposit has a recorded `sha256sum`
- [ ] Repository result footprint after cleanup measured; target 12–15 MB
- [ ] Retained in Git: `features/*.csv`, `scent_validation/` summaries and 3 PNGs,
      `scent_validation_min_distance/*.csv`, `scent_chr_sweep_summary.csv`, and per mode
      `summary_metrics`, `tier_summary`, `top100_links`, `distance_distribution.png`, and three
      `validation_*` summaries
- [ ] Duplicated baseline artifacts deduplicated — 6 patterns × 10–11 redundant copies, ~20.6 MB
      (`TODO.md` §4.1)
- [ ] Every archived path has a pointer in `docs/input_output_reference.md` §8
- [ ] `scent_validation_all_ranked_methods_combined.csv` resolved: Git LFS, or the min-distance
      rule made to depend on the validation rule rather than the committed file. It has a live
      downstream consumer, so removing it from the clone breaks re-running the controls
- [ ] `results/pbmc/audits/reranker_output_audit_checks.csv` reviewed — present in the tree,
      excluded from the archive, contents unknown
- [ ] A figure produced for the min-distance controls; none currently exists
      (`docs/results_report.md` §6)
- [ ] The two degenerate distance bins (`200_500kb`, `gt500kb`) either emit `NA` or are annotated
      as untested in the committed CSV, so the 8.906 artifact cannot be quoted from the raw file

---

## 10. What not to include in the public release

- [ ] `data/` — raw input, licensing unresolved
- [ ] `resources/jaspar/JASPAR2022.sqlite` — third-party, licensing unresolved
- [ ] `renv/library/`, `renv/staging/` — restored inside the image
- [ ] `tarballs/` — nested duplicate of `results/`, ~half the repository tree
- [ ] `results/pbmc/rankings_backup/` — superseded 7-mode copy predating the moddist work
- [ ] `results/pbmc.before_restore_20260715_195026/` — restore scratch
- [ ] All `*.bak` files, `run_analysis_before_scent.py`, `renv.lock.before_scent`
- [ ] `scent_sweep_file_inventory.txt`, `repo_tree_for_claude.txt`
- [ ] `config/config.yaml`, `containers/environment.yml` — both 0 bytes
- [ ] `README_SPLIT.md` — stale; describes the current architecture as a proposal
- [ ] `docs/*.docx` — three binary files, one with spaces in the filename. Extract, merge, delete
- [ ] Any absolute local paths, hostnames or usernames. Grep `scripts/` and `docs/` before
      publishing
- [ ] Any credentials or tokens in `docs/lab_notebook.md` — 66 KB of unedited notes, worth a
      `git grep` for `token`, `key`, `password`, `ssh`

### Decide explicitly

- [ ] **`docs/lab_notebook.md` (formerly `developer_notes.md`).** It contains candid negative
      self-assessment: "not publishable as-is", "Stop polishing this reranker",
      "Option A — Archive / stop", plus the ArchR integration failures. Recommendation: **publish
      it, with a header** stating it is an unedited chronological notebook and that
      `docs/results_report.md` is the current statement. It is the most honest artifact in the
      repository and it demonstrates the reasoning that produced the controls. Publishing it
      unlabelled, however, invites someone to quote a superseded verdict as the conclusion
- [ ] **`docs/benchmark_summary.md`.** Same issue, smaller. Keep, with a header noting it is a
      chronological record containing revised conclusions
- [ ] `scripts/dev/` — 15 exploratory scripts, not part of the workflow. Keep with a `README`
      stating they are unmaintained, or move to a branch

---

## 11. What should be archived externally

| Item | Size | Destination |
|---|---|---|
| 11 × `pbmc_*_ranked_links.csv` | ~29 MB | Zenodo — regenerable in minutes from the feature table |
| 22 × `scent_links_chr*.csv` | ~15 MB | Zenodo — expensive to regenerate; SCENT is the slow step |
| 22 × `scent_candidates_chr*.csv` | ~2 MB | Zenodo — **the 117,811-pair de novo candidate universe. Do not lose** |
| `scent_links_all_chromosomes.csv` | 6.7 MB | Zenodo |
| `scent_validation_all_ranked_methods_combined.csv` | 6.1 MB | Git LFS preferred — has a live downstream consumer |
| Per-mode `ora_dotplot.png`, `linkpeaks_vs_model_scatter.png` | ~7 MB | Zenodo |
| Per-mode `validation_promoted_demoted_inspection.csv`, promoted/demoted tables | ~2 MB | Zenodo |
| `resources/jaspar/JASPAR2022.sqlite` | large | LFS or runtime fetch, subject to §1 |
| `tarballs/*.tar.gz` | large | Zenodo if preserved at all, otherwise delete after §8 tagging |
| `data/` | large | Only if licensing permits; otherwise cite the accession |

---

## 12. Final pre-publication pass

- [ ] `README.md` first screen contains no claim stronger than the evidence. Specifically: no
      "outperforms" without the distance-only caveat, no "improves peak–gene linking", no
      "validated"
- [ ] `docs/results_report.md` §11 ("What cannot be claimed") is intact and linked from
      `README.md`
- [ ] `CITATION.cff` abstract says "benchmark", not "method"
- [ ] Repository name and description say "benchmark"
- [ ] The 8.906 odds-ratio artifact does not appear as a result anywhere in `docs/`
- [ ] Every internal documentation link resolves
- [ ] Every badge resolves
- [ ] Zenodo DOI resolves and the deposit matches `MANIFEST.sha256`
- [ ] Docker image pulls from a clean machine and the smoke test passes
- [ ] A clean clone plus `data/` and the JASPAR sqlite reproduces `--dry-run` for every section
- [ ] `docs/lab_notebook.md` and `docs/benchmark_summary.md` carry their headers
- [ ] All 15 items in `TODO.md` §9 are either resolved or explicitly documented as unknown in
      `README.md`
- [ ] `full` (λ = 0.3) described as an aggressive distance-prior sensitivity setting, never as
      the primary model or a headline result
- [ ] No claim that λ = 0.3 beats λ = 0.1 in a distance-controlled sense — the within-bin odds
      ratios are identical