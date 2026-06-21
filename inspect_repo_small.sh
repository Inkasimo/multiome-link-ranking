#!/usr/bin/env bash
set -euo pipefail

OUT="repo_inspection_small.txt"

{
  echo "============================================================"
  echo "BASIC"
  echo "============================================================"
  echo "PWD: $(pwd)"
  echo "BRANCH: $(git branch --show-current 2>/dev/null || echo NA)"
  echo
  echo "REMOTE:"
  git remote -v || true
  echo
  echo "STATUS:"
  git status --short --branch || true
  echo
  echo "RECENT COMMITS:"
  git log --oneline --decorate -n 5 || true
  echo

  echo "============================================================"
  echo "TOP LEVEL"
  echo "============================================================"
  find . -maxdepth 1 -mindepth 1 \
    -not -name ".git" \
    -printf "%f\n" | sort
  echo

  echo "============================================================"
  echo "IMPORTANT FILES"
  echo "============================================================"
  for f in \
    containers/Dockerfile \
    workflow/Snakefile \
    config/config.yaml \
    config/pbmc.yaml \
    run_analysis.py \
    wrapper-requirements.txt \
    .Rprofile \
    renv.lock
  do
    if [ -f "$f" ]; then
      echo "FOUND   $f"
    else
      echo "MISSING $f"
    fi
  done
  echo

  echo "============================================================"
  echo "TRACKED FILES"
  echo "============================================================"
  git ls-files | sort
  echo

  echo "============================================================"
  echo "SCRIPTS"
  echo "============================================================"
  if [ -d scripts ]; then
    find scripts -maxdepth 1 -type f -printf "%f\n" | sort
  else
    echo "MISSING scripts/"
  fi
  echo

  echo "============================================================"
  echo "CONTAINER / WORKFLOW PREVIEWS"
  echo "============================================================"
  for f in containers/Dockerfile workflow/Snakefile config/config.yaml run_analysis.py .Rprofile; do
    echo
    echo "----- $f -----"
    if [ -f "$f" ]; then
      sed -n '1,80p' "$f"
    else
      echo "MISSING"
    fi
  done

} > "$OUT"

echo
echo "Wrote full inspection to: $OUT"
echo
echo "==================== SCREEN SUMMARY ===================="
grep -E "^(PWD:|BRANCH:|FOUND|MISSING|##|On branch|STATUS:|RECENT COMMITS:)" "$OUT" || true
echo
echo "Tracked files:"
git ls-files | sort | sed -n '1,80p'
echo
echo "Scripts:"
if [ -d scripts ]; then
  find scripts -maxdepth 1 -type f -printf "%f\n" | sort
else
  echo "MISSING scripts/"
fi
echo
echo "Now paste the file with:"
echo "cat repo_inspection_small.txt"
