#!/usr/bin/env python3

import argparse
import subprocess
from pathlib import Path


def run(cmd):
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(
        description="Wrapper for the multiome link-ranking Snakemake workflow."
    )
    parser.add_argument("--cores", type=int, default=4)
    parser.add_argument("--snakefile", default="workflow/Snakefile")
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    cmd = [
        "snakemake",
        "-s", args.snakefile,
        "--configfile", args.config,
        "--cores", str(args.cores),
    ]

    if args.dry_run:
        cmd.append("-n")

    run(cmd)


if __name__ == "__main__":
    main()
