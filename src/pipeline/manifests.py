from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any, Iterable


def write_csv(path: str | Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k) for k in fieldnames})


def read_csv(path: str | Path) -> list[dict[str, Any]]:
    p = Path(path)
    if not p.exists():
        return []
    with p.open("r", newline="") as f:
        reader = csv.DictReader(f)
        return [dict(r) for r in reader]


def extract_design_id(name: str) -> int | None:
    tokens = re.split(r"[_-]+", name)
    if tokens and tokens[-1].isdigit():
        return int(tokens[-1])
    return None


def structure_id_from_name(name: str) -> str:
    tokens = re.split(r"[_-]+", name)
    if tokens and tokens[-1].isdigit():
        return "_".join(tokens[:-1]) or name
    return name


def find_metrics_file(output_dir: str | Path) -> Path | None:
    out = Path(output_dir)
    for cand in [
        out / "metrics.csv",
        out / "af3score_metrics.csv",
        out / "af3score_base_outputs" / "af3score_metrics.csv",
    ]:
        if cand.exists():
            return cand
    return None


def build_name_map(root: str | Path, pattern: str = "*.pdb") -> dict[str, Path]:
    mapping: dict[str, Path] = {}
    for fp in Path(root).rglob(pattern):
        mapping[fp.stem] = fp
    return mapping
