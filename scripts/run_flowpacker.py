#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import os
import shutil
import subprocess
import time
from pathlib import Path


def _parse_fasta(path: Path) -> list[str]:
    seqs: list[str] = []
    current: list[str] = []
    started = False
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if started:
                seqs.append("".join(current))
            current = []
            started = True
        else:
            current.append(line)
    if started:
        seqs.append("".join(current))
    return seqs


def build_seq_csv(seq_dir: Path, output_csv: Path, suffix: str) -> int:
    seen = set()
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with output_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["link_name", "seq", "seq_idx"])
        for fasta in sorted(seq_dir.rglob("*.fa*")):
            base_name = fasta.stem
            seqs = _parse_fasta(fasta)
            for i, seq in enumerate(seqs):
                if i == 0:
                    continue
                if not seq:
                    continue
                if seq in seen:
                    continue
                seen.add(seq)
                link_name = f"{base_name}{suffix}"
                writer.writerow([link_name, seq, str(i)])
                count += 1
    return count


def split_batches(pdbs: list[Path], out_dir: Path, num_jobs: int) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    if num_jobs <= 0:
        num_jobs = 1
    batch_size = math.ceil(len(pdbs) / num_jobs)
    batches: list[Path] = []
    for idx in range(num_jobs):
        batch = pdbs[idx * batch_size : (idx + 1) * batch_size]
        if not batch:
            continue
        batch_dir = out_dir / f"batch_{idx + 1}"
        batch_dir.mkdir(parents=True, exist_ok=True)
        for pdb in batch:
            dst = batch_dir / pdb.name
            if not dst.exists():
                try:
                    os.link(pdb, dst)
                except OSError:
                    shutil.copy2(pdb, dst)
        batches.append(batch_dir)
    return batches


def _abspath_if_relative(value: str | None, root: Path) -> str | None:
    if not value:
        return value
    p = Path(value)
    if p.is_absolute():
        return str(p)
    return str((root / p).resolve())


def write_batch_yamls(base_yaml: Path, batch_dirs: list[Path], out_dir: Path, flowpacker_repo: Path) -> list[Path]:
    import yaml

    out_dir.mkdir(parents=True, exist_ok=True)
    base_cfg = yaml.safe_load(base_yaml.read_text())
    # Make ckpt paths absolute so downstream sampler works from any cwd.
    if isinstance(base_cfg, dict):
        base_cfg["ckpt"] = _abspath_if_relative(base_cfg.get("ckpt"), flowpacker_repo)
        base_cfg["conf_ckpt"] = _abspath_if_relative(base_cfg.get("conf_ckpt"), flowpacker_repo)
    yaml_paths: list[Path] = []
    for batch_dir in batch_dirs:
        cfg = dict(base_cfg)
        cfg.setdefault("data", {})
        cfg["data"]["test_path"] = str(batch_dir)
        yaml_path = out_dir / f"{batch_dir.name}.yml"
        yaml_path.write_text(yaml.safe_dump(cfg, sort_keys=False))
        yaml_paths.append(yaml_path)
    return yaml_paths


def _resolve_base_yaml(flowpacker_repo: Path, base_yaml: str | None) -> Path:
    if base_yaml:
        return Path(base_yaml)
    candidates = [
        flowpacker_repo / "config" / "inference" / "base.yaml",
        flowpacker_repo / "config" / "inference" / "base.yml",
        flowpacker_repo / "config" / "base.yaml",
        flowpacker_repo / "config" / "base.yml",
    ]
    for cand in candidates:
        if cand.exists():
            return cand
    raise FileNotFoundError("Could not locate FlowPacker base.yaml; pass --base_yaml explicitly.")


def _progress(idx: int, total: int, item: str, status: str, elapsed: float) -> None:
    print(
        f"[flowpacker] {idx}/{total} {status} elapsed={elapsed:.2f}s item={item}",
        flush=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_pdb_dir", required=True)
    parser.add_argument("--seq_fasta_dir", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--flowpacker_repo", required=True)
    parser.add_argument("--base_yaml", default=None)
    parser.add_argument("--num_jobs", type=int, default=1)
    parser.add_argument("--binder_chain", default="A")
    parser.add_argument("--link_suffix", default=".pdb")
    args = parser.parse_args()

    input_pdb_dir = Path(args.input_pdb_dir).resolve()
    seq_fasta_dir = Path(args.seq_fasta_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    flowpacker_repo = Path(args.flowpacker_repo).resolve()

    if not input_pdb_dir.exists():
        raise FileNotFoundError(f"input_pdb_dir not found: {input_pdb_dir}")
    if not seq_fasta_dir.exists():
        raise FileNotFoundError(f"seq_fasta_dir not found: {seq_fasta_dir}")
    if not flowpacker_repo.exists():
        raise FileNotFoundError(f"flowpacker_repo not found: {flowpacker_repo}")

    csv_path = output_dir / "flowpacker_input.csv"
    seq_count = build_seq_csv(seq_fasta_dir, csv_path, args.link_suffix)
    if seq_count == 0:
        raise RuntimeError(f"No sequences written to {csv_path}")

    pdbs = sorted(input_pdb_dir.glob("*.pdb"))
    if not pdbs:
        pdbs = sorted(input_pdb_dir.rglob("*.pdb"))
    if not pdbs:
        raise RuntimeError(f"No PDBs found in {input_pdb_dir}")

    batch_dir = output_dir / "input_pdb_batch"
    batches = split_batches(pdbs, batch_dir, args.num_jobs)

    base_yaml = _resolve_base_yaml(flowpacker_repo, args.base_yaml)
    yaml_dir = output_dir / "batch_yml"
    yaml_paths = write_batch_yamls(base_yaml, batches, yaml_dir, flowpacker_repo)

    flowpacker_out = output_dir / "flowpacker_outputs"
    flowpacker_out.mkdir(parents=True, exist_ok=True)

    sampler_script = Path(__file__).resolve().parent / "flowpacker_sampler_pipe.py"
    if not sampler_script.exists():
        raise FileNotFoundError(f"Missing sampler script: {sampler_script}")

    env = os.environ.copy()
    env["FLOWPACKER_REPO"] = str(flowpacker_repo)
    env["FLOWPACKER_BINDER_CHAIN"] = args.binder_chain

    total = len(yaml_paths)
    for idx, yaml_path in enumerate(yaml_paths, start=1):
        cmd = [
            "python",
            str(sampler_script),
            str(yaml_path),
            "--save_dir",
            str(flowpacker_out),
            "--csv_file",
            str(csv_path),
            "--binder_chain",
            args.binder_chain,
        ]
        start = time.time()
        status = "OK"
        try:
            subprocess.check_call(cmd, env=env)
        except Exception:
            status = "FAILED"
            raise
        finally:
            _progress(idx, total, yaml_path.stem, status, time.time() - start)


if __name__ == "__main__":
    main()
