#!/usr/bin/env bash
#
# download_inputs.sh — prepare all inputs required by the workflow.
#
# Fetches the 10x Genomics pbmc_unsorted_10k multiome files into data/,
# stripping the sample-name prefix so the filenames match config/default.yaml.
# Verifies every file against resources/input_manifest.tsv.
#
# JASPAR2022.sqlite is NOT fetched from upstream: the jaspar2022.genereg.net
# URL referenced in scripts/run_linkpeaks_reranker.R is a BiocFileCache key,
# not a live download target, and JASPAR has since moved to ELIXIR hosting.
# The file is distributed with the release archive instead.
#
# Usage:
#   bash scripts/download_inputs.sh            # fetch missing, verify all
#   bash scripts/download_inputs.sh --verify   # verify only, never download
#
# Run from the repository root.

set -euo pipefail

MANIFEST="resources/input_manifest.tsv"
VERIFY_ONLY=0
[ "${1:-}" = "--verify" ] && VERIFY_ONLY=1

# Fill this in after the Zenodo deposit is published, then commit.
JASPAR_SOURCE_URL="https://zenodo.org/records/22032568/files/JASPAR2022.sqlite?download=1"

if [ ! -f "workflow/Snakefile" ] || [ ! -f "config/default.yaml" ]; then
  echo "ERROR: run this from the repository root." >&2
  exit 2
fi

if [ ! -f "$MANIFEST" ]; then
  echo "ERROR: manifest not found: $MANIFEST" >&2
  exit 2
fi

command -v sha256sum >/dev/null || { echo "ERROR: sha256sum not found." >&2; exit 2; }

fetch() {
  # fetch <url> <destination>
  if command -v curl >/dev/null; then
    curl -fL --retry 3 --continue-at - -o "$2" "$1"
  elif command -v wget >/dev/null; then
    wget -c -O "$2" "$1"
  else
    echo "ERROR: neither curl nor wget available." >&2
    exit 2
  fi
}

verify() {
  # verify <path> <expected_sha256>  -> 0 if match
  echo "$2  $1" | sha256sum --check --status 2>/dev/null
}

rc=0

while IFS=$'\t' read -r path sha bytes source; do
  case "$path" in
    ''|'#'*) continue ;;
  esac

  # ---- already present and correct -------------------------------------
  if [ -f "$path" ]; then
    if verify "$path" "$sha"; then
      echo "OK       $path"
      continue
    fi
    echo "MISMATCH $path" >&2
    echo "         expected $sha" >&2
    echo "         got      $(sha256sum "$path" | cut -d' ' -f1)" >&2
    echo "         delete the file and re-run to re-fetch." >&2
    rc=1
    continue
  fi

  # ---- missing ----------------------------------------------------------
  if [ "$VERIFY_ONLY" -eq 1 ]; then
    echo "MISSING  $path" >&2
    rc=1
    continue
  fi

  if [ "$source" = "ZENODO_PENDING" ]; then
    echo "MISSING  $path" >&2
    if [ -n "$JASPAR_SOURCE_URL" ]; then
      echo "         fetching from the release archive..." >&2
      mkdir -p "$(dirname "$path")"
      fetch "$JASPAR_SOURCE_URL" "$path"
      if verify "$path" "$sha"; then
        echo "OK       $path (downloaded)"
        continue
      fi
      echo "ERROR: checksum mismatch after download: $path" >&2
      rc=1
      continue
    fi
    echo "         JASPAR2022.sqlite ships with the release archive, not with Git." >&2
    echo "         Obtain it from the Zenodo deposit linked in README.md, place it at" >&2
    echo "         $path, then re-run this script." >&2
    echo "         Expected sha256: $sha" >&2
    rc=1
    continue
  fi

  echo "FETCH    $path" >&2
  echo "         <- $source" >&2
  mkdir -p "$(dirname "$path")"
  fetch "$source" "$path"

  if verify "$path" "$sha"; then
    echo "OK       $path (downloaded)"
  else
    echo "ERROR: checksum mismatch after download: $path" >&2
    echo "       expected $sha" >&2
    echo "       got      $(sha256sum "$path" | cut -d' ' -f1)" >&2
    rc=1
  fi
done < "$MANIFEST"

echo
if [ "$rc" -eq 0 ]; then
  echo "All inputs present and verified."
else
  echo "One or more inputs are missing or failed verification. See messages above." >&2
fi

exit "$rc"
