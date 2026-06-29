#!/usr/bin/env python3

import argparse
import os
import shlex
import subprocess
import sys
from pathlib import Path

import yaml


DEFAULT_IMAGE = "multiome-link-ranking:pilot"


def q(x: str) -> str:
    return shlex.quote(str(x))


def print_cmd(cmd: list[str]) -> None:
    print("Running:\n  " + " \\\n  ".join(q(c) for c in cmd), flush=True)


def run(cmd: list[str]) -> int:
    print_cmd(cmd)
    return subprocess.call(cmd)


def load_yaml(path: str) -> dict:
    with open(path) as handle:
        out = yaml.safe_load(handle)
    return out or {}


def docker_bind_mount_works(repo_root: Path, image: str) -> bool:
    cmd = [
        "docker", "run", "--rm",
        "-v", f"{repo_root}:/work",
        "-w", "/work",
        image,
        "sh", "-lc",
        "test -f workflow/Snakefile && test -f config/default.yaml",
    ]

    return subprocess.call(
        cmd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    ) == 0


def build_common_snakemake_args(args) -> list[str]:
    smk = [
        "snakemake",
        "-s", args.snakefile,
        "--configfile", args.configfile,
        "--cores", str(args.cores),
        "--reason",
        "-p",
    ]

    if args.rerun_incomplete:
        smk.append("--rerun-incomplete")

    if args.rerun_triggers:
        smk.extend(["--rerun-triggers", args.rerun_triggers])

    if args.dry_run:
        smk.append("-n")

    if args.extra:
        extra = args.extra
        if extra and extra[0] == "--":
            extra = extra[1:]
        smk.extend(extra)

    return smk


def main() -> int:
    parser = argparse.ArgumentParser(
        prog="run_analysis.py",
        description="Docker wrapper for the LinkPeaks reranker scoring workflow.",
    )

    parser.add_argument(
        "section",
        choices=[
            "build_linkpeaks_features",
            "run_default_score",
            "run_score_mode",
            "run_all_score_modes",
            "run_reranker_score_suite",
            "run_scent_sweep",
            "run_scent_validation",
            "run_scent_pipeline",
            "run_reranker_with_scent",
            "list_score_modes",
            "list_scent_run",
            "list_scent_methods",
            "unlock",
        ],
        help=(
            "What to run: build_linkpeaks_features, run_default_score, "
            "run_score_mode, run_all_score_modes, run_reranker_score_suite, "
            "run_scent_sweep, run_scent_validation, run_scent_pipeline, "
            "run_reranker_with_scent, list_score_modes, list_scent_run, "
            "list_scent_methods, unlock."
        ),
    )

    parser.add_argument(
        "--mode",
        default=None,
        help=(
            "Score mode for run_score_mode, for example: "
            "full, coactivity, distance_only, full_linkpeaks_anchored."
        ),
    )

    parser.add_argument("--image", default=DEFAULT_IMAGE)
    parser.add_argument("--snakefile", default="workflow/Snakefile")
    parser.add_argument("--configfile", default="config/default.yaml")

    parser.add_argument("--cores", type=int, default=4)
    parser.add_argument("--cpus", type=int, default=4)

    parser.add_argument("--dry-run", action="store_true")

    parser.add_argument(
        "--rerun-incomplete",
        dest="rerun_incomplete",
        action="store_true",
        default=True,
        help="Rerun incomplete Snakemake jobs. Default: on.",
    )
    parser.add_argument(
        "--no-rerun-incomplete",
        dest="rerun_incomplete",
        action="store_false",
        help="Disable Snakemake --rerun-incomplete.",
    )

    parser.add_argument(
        "--rerun-triggers",
        default="mtime",
        help="Snakemake rerun trigger policy. Default: mtime.",
    )

    parser.add_argument(
        "--extra",
        nargs=argparse.REMAINDER,
        help="Extra arguments passed to Snakemake. Use this last.",
    )

    args = parser.parse_args()
    repo_root = Path.cwd()

    if not (repo_root / args.snakefile).exists():
        print(f"ERROR: Snakefile not found: {args.snakefile}", file=sys.stderr)
        return 2

    if not (repo_root / args.configfile).exists():
        print(f"ERROR: config file not found: {args.configfile}", file=sys.stderr)
        return 2

    try:
        cfg = load_yaml(args.configfile)
    except Exception as e:
        print(f"ERROR: could not read config file {args.configfile}: {e}", file=sys.stderr)
        return 2

    dataset = cfg["dataset"]
    output_root = cfg.get("output_root", "results")

    score_modes_file = cfg.get(
        "score_modes_file",
        cfg.get("ablations_file", "config/ablations.yaml"),
    )

    default_score_mode = cfg.get("default_score_mode", "full")

    if not (repo_root / score_modes_file).exists():
        print(f"ERROR: score modes file not found: {score_modes_file}", file=sys.stderr)
        return 2

    try:
        score_modes_cfg = load_yaml(score_modes_file)
    except Exception as e:
        print(f"ERROR: could not read score modes file {score_modes_file}: {e}", file=sys.stderr)
        return 2

    score_modes = score_modes_cfg.get("score_modes")
    if score_modes is None:
        score_modes = score_modes_cfg.get("ablations", {})

    modes = list(score_modes.keys())

    if not modes:
        print(f"ERROR: no score modes found in {score_modes_file}", file=sys.stderr)
        return 2

    if default_score_mode not in modes:
        print(
            f"ERROR: default_score_mode='{default_score_mode}' is not present in {score_modes_file}",
            file=sys.stderr,
        )
        print("Available modes: " + ", ".join(modes), file=sys.stderr)
        return 2

    features_done = f"{output_root}/{dataset}/features/.done"

    all_score_mode_targets = [
        f"{output_root}/{dataset}/rankings/{mode}/.done"
        for mode in modes
    ]

    scent_run_config_file = cfg.get("scent_run_config", "config/scent_run.yaml")
    scent_run_cfg = {}
    if (repo_root / scent_run_config_file).exists():
        try:
            scent_run_raw = load_yaml(scent_run_config_file)
            scent_run_cfg = scent_run_raw.get("scent_run", scent_run_raw)
        except Exception as e:
            print(f"ERROR: could not read SCENT run config {scent_run_config_file}: {e}", file=sys.stderr)
            return 2

    scent_sweep_dir = scent_run_cfg.get(
        "output_dir",
        f"{output_root}/{dataset}/scent_chr_sweep_100kb_frac020_1000cells",
    )
    scent_sweep_done = f"{scent_sweep_dir}/.done"

    scent_config_file = cfg.get("scent_validation_config", "config/scent_validation.yaml")
    scent_cfg = {}
    if (repo_root / scent_config_file).exists():
        try:
            scent_raw = load_yaml(scent_config_file)
            scent_cfg = scent_raw.get("scent_validation", scent_raw)
        except Exception as e:
            print(f"ERROR: could not read SCENT validation config {scent_config_file}: {e}", file=sys.stderr)
            return 2

    scent_methods = [str(x) for x in scent_cfg.get(
        "methods",
        ["linkpeaks", "coactivity", "coactivity_tf", "full_moddist_lambda_0_1", "distance_only"],
    )]
    scent_output_dir = scent_cfg.get(
        "output_dir",
        f"{output_root}/{dataset}/scent_validation",
    )
    scent_done = f"{scent_output_dir}/.done"

    if args.section == "list_score_modes":
        print(f"Default score mode: {default_score_mode}")
        print("Available score modes:")
        for mode in modes:
            entry = score_modes[mode] or {}
            score_mode = entry.get("score_mode", mode)
            lambda_distance = entry.get("lambda_distance", cfg.get("lambda_distance"))
            alpha_tf = entry.get("alpha_tf", cfg.get("alpha_tf"))

            label = " [default]" if mode == default_score_mode else ""
            print(
                f"  - {mode}{label} "
                f"(score_mode={score_mode}, "
                f"lambda_distance={lambda_distance}, "
                f"alpha_tf={alpha_tf})"
            )
        return 0

    if args.section == "list_scent_run":
        print(f"SCENT run config: {scent_run_config_file}")
        print(f"SCENT sweep done target: {scent_sweep_done}")
        print(f"SCENT sweep output directory: {scent_sweep_dir}")
        print(f"Chromosomes: {scent_run_cfg.get('chromosomes', ['chr1'])}")
        print(f"link_distance: {scent_run_cfg.get('link_distance', cfg.get('link_distance', 100000))}")
        print(f"min_pair_frac: {scent_run_cfg.get('min_pair_frac', 0.02)}")
        print(f"max_cells: {scent_run_cfg.get('max_cells', 1000)}")
        print(f"max_scent_candidates: {scent_run_cfg.get('max_scent_candidates', 100000)}")
        print(f"scent_cores: {scent_run_cfg.get('scent_cores', 4)}")
        print(f"scent_regr: {scent_run_cfg.get('scent_regr', 'poisson')}")
        return 0

    if args.section == "list_scent_methods":
        print(f"SCENT validation config: {scent_config_file}")
        print(f"SCENT output done target: {scent_done}")
        print("SCENT validation methods:")
        missing = []
        for mode in scent_methods:
            status = "" if mode in modes else " [MISSING from score modes]"
            if mode not in modes:
                missing.append(mode)
            print(f"  - {mode}{status}")
        if missing:
            print("ERROR: these SCENT methods are not present in the ablations/score modes file:", file=sys.stderr)
            print("  " + ", ".join(missing), file=sys.stderr)
            return 2
        return 0

    if not docker_bind_mount_works(repo_root, args.image):
        print(
            "ERROR: Docker bind mount failed.\n\n"
            f"Host path: {repo_root}\n"
            "The container could not see /work/workflow/Snakefile.\n"
            "Check Docker Desktop file-sharing / WSL mount access.",
            file=sys.stderr,
        )
        return 2

    docker = [
        "docker", "run", "--rm", "-i",
        "--user", f"{os.getuid()}:{os.getgid()}",
        "-e", "HOME=/tmp",
        "-e", "XDG_CACHE_HOME=/tmp/.cache",
        "-e", "XDG_CONFIG_HOME=/tmp/.config",
        "--init",
        "--cpus", str(args.cpus),
        "-v", f"{repo_root}:/work",
        "-w", "/work",
        args.image,
    ]

    if args.section == "unlock":
        smk = [
            "snakemake",
            "-s", args.snakefile,
            "--configfile", args.configfile,
            "--unlock",
        ]
        return run(docker + smk)

    targets: list[str] = []

    if args.section == "build_linkpeaks_features":
        targets = [features_done]

    elif args.section == "run_default_score":
        targets = [f"{output_root}/{dataset}/rankings/{default_score_mode}/.done"]

    elif args.section == "run_score_mode":
        if args.mode is None:
            print("ERROR: run_score_mode requires --mode.", file=sys.stderr)
            print("Available modes: " + ", ".join(modes), file=sys.stderr)
            return 2

        if args.mode not in modes:
            print(f"ERROR: unknown score mode: {args.mode}", file=sys.stderr)
            print("Available modes: " + ", ".join(modes), file=sys.stderr)
            return 2

        targets = [f"{output_root}/{dataset}/rankings/{args.mode}/.done"]

    elif args.section == "run_all_score_modes":
        targets = all_score_mode_targets

    elif args.section == "run_reranker_score_suite":
        targets = ["all"]

    elif args.section == "run_scent_sweep":
        targets = [scent_sweep_done]

    elif args.section == "run_scent_validation":
        missing = [m for m in scent_methods if m not in modes]
        if missing:
            print("ERROR: SCENT validation methods missing from score modes: " + ", ".join(missing), file=sys.stderr)
            return 2
        targets = [scent_done]

    elif args.section == "run_scent_pipeline":
        missing = [m for m in scent_methods if m not in modes]
        if missing:
            print("ERROR: SCENT validation methods missing from score modes: " + ", ".join(missing), file=sys.stderr)
            return 2
        targets = [scent_sweep_done, scent_done]

    elif args.section == "run_reranker_with_scent":
        missing = [m for m in scent_methods if m not in modes]
        if missing:
            print("ERROR: SCENT validation methods missing from score modes: " + ", ".join(missing), file=sys.stderr)
            return 2
        targets = ["all", scent_done]

    else:
        print(f"ERROR: unhandled section: {args.section}", file=sys.stderr)
        return 2

    smk = build_common_snakemake_args(args)

    if targets:
        smk.append("--")
        smk.extend(targets)

    return run(docker + smk)


if __name__ == "__main__":
    raise SystemExit(main())
