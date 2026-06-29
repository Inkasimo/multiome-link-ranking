#!/usr/bin/env bash
set -euo pipefail

DATASET="${1:-pbmc}"
OUTDIR="tarballs"
STAMP="$(date +%Y%m%d_%H%M%S)"
TARBALL="${OUTDIR}/${DATASET}_reranker_validation_${STAMP}.tar.gz"
MANIFEST="${OUTDIR}/${DATASET}_reranker_validation_${STAMP}_manifest.txt"

mkdir -p "$OUTDIR"

echo "Creating validation tarball for dataset: ${DATASET}"
echo "Output: ${TARBALL}"

{
  echo "Reranker validation export"
  echo "Dataset: ${DATASET}"
  echo "Created: $(date -Is)"
  echo
  echo "Included purpose:"
  echo "- reranking outputs"
  echo "- validation summaries"
  echo "- distance/diversity diagnostics"
  echo "- modified-distance ablation outputs"
  echo "- configs and scoring code needed to inspect behavior"
  echo
  echo "Score mode directories present:"
  if [ -d "results/${DATASET}/rankings" ]; then
    find "results/${DATASET}/rankings" -maxdepth 1 -mindepth 1 -type d | sort
  else
    echo "MISSING: results/${DATASET}/rankings"
  fi
  echo
  echo "Feature files present:"
  find "results/${DATASET}/features" -maxdepth 1 -type f 2>/dev/null | sort || true
  echo
  echo "Key validation files:"
  find "results/${DATASET}/rankings" -type f \
    \( -name "*summary_metrics.csv" \
       -o -name "*validation*.csv" \
       -o -name "*topN_overlap*.csv" \
       -o -name "*top100_links.csv" \
       -o -name "*top_promoted*.csv" \
       -o -name "*top_demoted*.csv" \
       -o -name "*ranked_links.csv" \
       -o -name "*.png" \) \
    2>/dev/null | sort || true
} > "$MANIFEST"

tar -czf "$TARBALL" \
  "$MANIFEST" \
  config/default.yaml \
  config/ablations.yaml \
  workflow/Snakefile \
  run_analysis.py \
  Dockerfile \
  README.md \
  README_SPLIT.md \
  scripts/evaluate_rankings.R \
  scripts/run_linkpeaks_reranker.R \
  results/"${DATASET}"/features/*.csv \
  results/"${DATASET}"/rankings \
  2>/tmp/reranker_tar_warnings.txt || {
    echo "Tar failed. Warnings/errors:"
    cat /tmp/reranker_tar_warnings.txt
    exit 1
  }

echo
echo "Done:"
ls -lh "$TARBALL"
echo
echo "Manifest:"
ls -lh "$MANIFEST"
echo
echo "Tarball contents preview:"
tar -tzf "$TARBALL" | head -80
