#!/usr/bin/env bash
set -euo pipefail

echo "============================================================"
echo "REPO INSPECTION"
echo "============================================================"
echo

echo "PWD:"
pwd
echo

echo "DATE:"
date
echo

echo "GIT REMOTE:"
git remote -v || true
echo

echo "BRANCH:"
git branch --show-current || true
echo

echo "GIT STATUS:"
git status --short --branch || true
echo

echo "RECENT COMMITS:"
git log --oneline --decorate -n 10 || true
echo

echo "============================================================"
echo "TOP-LEVEL TREE"
echo "============================================================"
find . -maxdepth 2 \
  -path './.git' -prune -o \
  -path './results' -prune -o \
  -path './data' -prune -o \
  -print | sort
echo

echo "============================================================"
echo "TRACKED FILES"
echo "============================================================"
git ls-files | sort
echo

echo "============================================================"
echo "IGNORED FILES / DIRS"
echo "============================================================"
git status --ignored --short | sed -n '1,200p'
echo

echo "============================================================"
echo "KEY DIRECTORIES"
echo "============================================================"
for d in config workflow scripts containers docs renv results data; do
  echo
  echo "--- $d ---"
  if [ -d "$d" ]; then
    find "$d" -maxdepth 2 -type f | sort | sed -n '1,120p'
  else
    echo "MISSING"
  fi
done
echo

echo "============================================================"
echo "SCRIPT HEADERS"
echo "============================================================"
if [ -d scripts ]; then
  for f in scripts/*; do
    [ -f "$f" ] || continue
    echo
    echo "--- $f ---"
    sed -n '1,40p' "$f" || true
  done
fi
echo

echo "============================================================"
echo "CONFIG / WORKFLOW / CONTAINER PREVIEWS"
echo "============================================================"
for f in config/config.yaml config/pbmc.yaml workflow/Snakefile containers/Dockerfile wrapper-requirements.txt run_analysis.py .Rprofile renv.lock; do
  echo
  echo "--- $f ---"
  if [ -f "$f" ]; then
    sed -n '1,120p' "$f" || true
  else
    echo "MISSING"
  fi
done
echo

echo "============================================================"
echo "RESULTS OVERVIEW"
echo "============================================================"
if [ -d results ]; then
  find results -maxdepth 3 -type f | sort | sed -n '1,200p'
else
  echo "No results directory"
fi
echo

echo "============================================================"
echo "SCENT SWEEP OVERVIEW"
echo "============================================================"
if [ -d results/scent_chr_sweep_100kb_frac020_1000cells ]; then
  find results/scent_chr_sweep_100kb_frac020_1000cells -maxdepth 2 -type f | sort
else
  echo "No SCENT sweep directory"
fi
echo

echo "============================================================"
echo "BENCHMARK OUTPUT OVERVIEW"
echo "============================================================"
if [ -d results/benchmark_reranker_scent_sweep ]; then
  find results/benchmark_reranker_scent_sweep -maxdepth 1 -type f | sort
else
  echo "No benchmark output directory"
fi
echo

echo "============================================================"
echo "DONE"
echo "============================================================"
