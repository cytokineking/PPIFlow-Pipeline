#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def _dockq_cmd(dockq_bin: str) -> list[str]:
    if dockq_bin.endswith(".py"):
        return [sys.executable, dockq_bin]
    return [dockq_bin]


def main() -> int:
    parser = argparse.ArgumentParser(description="Run DockQ over a directory of PDBs")
    parser.add_argument("--dockq_bin", required=True, help="Path to DockQ executable or DockQ.py")
    parser.add_argument("--input_pdb_dir", required=True, help="Directory with model PDBs (e.g., relax outputs)")
    parser.add_argument("--reference_pdb_dir", required=True, help="Directory with reference PDBs")
    parser.add_argument("--output_dir", required=True, help="Directory to write *_dockq_score files")
    parser.add_argument("--allowed_mismatches", type=int, default=10, help="DockQ allowed mismatches")
    parser.add_argument("--skip_existing", action="store_true", help="Skip if output exists")

    args = parser.parse_args()
    input_dir = Path(args.input_pdb_dir)
    ref_dir = Path(args.reference_pdb_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    dockq_cmd = _dockq_cmd(args.dockq_bin)

    failures: list[str] = []
    for model_path in sorted(input_dir.glob("*.pdb")):
        ref_path = ref_dir / model_path.name
        if not ref_path.exists():
            failures.append(f"missing_reference\t{model_path.name}")
            continue
        out_path = out_dir / f"{model_path.stem}_dockq_score"
        if args.skip_existing and out_path.exists():
            continue
        cmd = (
            dockq_cmd
            + [
                "--allowed_mismatches",
                str(args.allowed_mismatches),
                str(model_path),
                str(ref_path),
                "--short",
            ]
        )
        try:
            with out_path.open("w") as handle:
                res = subprocess.run(cmd, stdout=handle, stderr=subprocess.PIPE, text=True)
            if res.returncode != 0:
                failures.append(f"dockq_failed\t{model_path.name}\t{res.stderr.strip()}")
        except Exception as exc:
            failures.append(f"dockq_exception\t{model_path.name}\t{exc}")

    if failures:
        fail_path = out_dir / "failed_records.txt"
        fail_path.write_text("\n".join(failures) + "\n")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
