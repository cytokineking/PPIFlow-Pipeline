#!/usr/bin/env python
from __future__ import annotations

import argparse
import ast
import csv
import math
import os
import subprocess
from pathlib import Path
from typing import Dict, Tuple


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _default_paths() -> Tuple[Path, Path]:
    root = _repo_root()
    og = root / "PPIFlow_OG_repo" / "demo_scripts" / "interface_analysis" / "get_interface_energy.py"
    new = root / "src" / "pipeline" / "rosetta" / "get_interface_energy.py"
    return og, new


def _run_parser(script: Path, input_pdbdir: Path, rosetta_dir: Path, binder_id: str, target_id: str, output_dir: Path, interface_dist: float, python_bin: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        python_bin,
        str(script),
        "--input_pdbdir",
        str(input_pdbdir),
        "--rosetta_dir",
        str(rosetta_dir),
        "--binder_id",
        binder_id,
        "--target_id",
        target_id,
        "--output_dir",
        str(output_dir),
        "--interface_dist",
        str(interface_dist),
    ]
    env = os.environ.copy()
    env.setdefault("MPLBACKEND", "Agg")
    subprocess.check_call(cmd, env=env)


def _load_binder_energy(csv_path: Path) -> Dict[str, Dict[str, float]]:
    results: Dict[str, Dict[str, float]] = {}
    with csv_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            name = row.get("pdbname") or Path(row.get("pdbpath", "")).stem
            raw = row.get("binder_energy") or "{}"
            try:
                parsed = ast.literal_eval(raw)
            except Exception:
                parsed = {}
            if not isinstance(parsed, dict):
                parsed = {}
            norm = {str(k): float(v) for k, v in parsed.items()}
            results[str(name)] = norm
    return results


def _summarize(d: Dict[str, Dict[str, float]]) -> Dict[str, float]:
    sums: Dict[str, float] = {}
    for k, v in d.items():
        sums[k] = sum(v.values())
    return sums


def _write_summary(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    og_default, new_default = _default_paths()

    parser = argparse.ArgumentParser(description="Compare OG vs new Rosetta interface parsing on identical outputs")
    parser.add_argument("--input_pdbdir", required=True, type=str, help="Directory of input PDBs")
    parser.add_argument("--rosetta_dir", required=True, type=str, help="Directory with Rosetta outputs per PDB")
    parser.add_argument("--binder_id", required=True, type=str, help="Binder chain id (e.g. A)")
    parser.add_argument("--target_id", required=True, type=str, help="Target chain id(s), comma-separated")
    parser.add_argument("--output_dir", required=True, type=str, help="Base output dir for comparison")
    parser.add_argument("--interface_dist", type=float, default=12.0, help="Interface distance cutoff")
    parser.add_argument("--og_script", type=str, default=str(og_default), help="Path to OG get_interface_energy.py")
    parser.add_argument("--new_script", type=str, default=str(new_default), help="Path to new get_interface_energy.py")
    parser.add_argument("--python", type=str, default=os.environ.get("PYTHON", "python"), help="Python executable")
    parser.add_argument("--tolerance", type=float, default=1e-6, help="Value diff tolerance")
    parser.add_argument("--skip_run", action="store_true", help="Skip running parsers and only compare existing CSVs")
    args = parser.parse_args()

    input_pdbdir = Path(args.input_pdbdir).resolve()
    rosetta_dir = Path(args.rosetta_dir).resolve()
    out_base = Path(args.output_dir).resolve()
    og_script = Path(args.og_script).resolve()
    new_script = Path(args.new_script).resolve()

    if not input_pdbdir.exists():
        raise SystemExit(f"input_pdbdir not found: {input_pdbdir}")
    if not rosetta_dir.exists():
        raise SystemExit(f"rosetta_dir not found: {rosetta_dir}")
    if not og_script.exists():
        raise SystemExit(f"OG script not found: {og_script}")
    if not new_script.exists():
        raise SystemExit(f"new script not found: {new_script}")

    og_out = out_base / "og"
    new_out = out_base / "new"

    if not args.skip_run:
        _run_parser(og_script, input_pdbdir, rosetta_dir, args.binder_id, args.target_id, og_out, args.interface_dist, args.python)
        _run_parser(new_script, input_pdbdir, rosetta_dir, args.binder_id, args.target_id, new_out, args.interface_dist, args.python)

    og_csv = og_out / "residue_energy.csv"
    new_csv = new_out / "residue_energy.csv"
    if not og_csv.exists():
        raise SystemExit(f"Missing OG residue_energy.csv at {og_csv}")
    if not new_csv.exists():
        raise SystemExit(f"Missing new residue_energy.csv at {new_csv}")

    og = _load_binder_energy(og_csv)
    new = _load_binder_energy(new_csv)
    og_sum = _summarize(og)
    new_sum = _summarize(new)

    all_names = sorted(set(og.keys()) | set(new.keys()))
    rows = []
    missing_og = 0
    missing_new = 0
    sum_diffs = []
    max_sum_diff = 0.0
    max_sum_name = ""

    for name in all_names:
        d1 = og.get(name)
        d2 = new.get(name)
        if d1 is None:
            missing_og += 1
            continue
        if d2 is None:
            missing_new += 1
            continue
        keys1 = set(d1.keys())
        keys2 = set(d2.keys())
        overlap = keys1 & keys2
        extra_og = keys1 - keys2
        extra_new = keys2 - keys1

        per_key_diffs = []
        for k in overlap:
            diff = abs(d1[k] - d2[k])
            if diff > args.tolerance:
                per_key_diffs.append(diff)
        max_key_diff = max(per_key_diffs) if per_key_diffs else 0.0
        mean_key_diff = sum(per_key_diffs) / len(per_key_diffs) if per_key_diffs else 0.0

        s1 = og_sum.get(name, 0.0)
        s2 = new_sum.get(name, 0.0)
        sum_diff = abs(s1 - s2)
        sum_diffs.append(sum_diff)
        if sum_diff > max_sum_diff:
            max_sum_diff = sum_diff
            max_sum_name = name

        rows.append({
            "pdb_name": name,
            "og_keys": len(keys1),
            "new_keys": len(keys2),
            "overlap_keys": len(overlap),
            "extra_og_keys": len(extra_og),
            "extra_new_keys": len(extra_new),
            "sum_og": s1,
            "sum_new": s2,
            "sum_abs_diff": sum_diff,
            "max_key_diff": max_key_diff,
            "mean_key_diff": mean_key_diff,
        })

    summary_path = out_base / "compare_summary.csv"
    _write_summary(summary_path, rows)

    def _stats(values: list[float]) -> Tuple[float, float, float]:
        if not values:
            return (0.0, 0.0, 0.0)
        v_sorted = sorted(values)
        mean = sum(values) / len(values)
        median = v_sorted[len(v_sorted) // 2]
        return (mean, median, max(values))

    mean_diff, median_diff, max_diff = _stats(sum_diffs)

    print("=== Rosetta interface parser comparison ===")
    print(f"Total PDBs: {len(all_names)}")
    print(f"Missing in OG: {missing_og}")
    print(f"Missing in new: {missing_new}")
    print(f"Sum diff mean: {mean_diff:.6f}")
    print(f"Sum diff median: {median_diff:.6f}")
    print(f"Sum diff max: {max_diff:.6f} ({max_sum_name})")
    print(f"Wrote per-PDB summary: {summary_path}")


if __name__ == "__main__":
    main()
