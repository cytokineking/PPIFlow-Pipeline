from __future__ import annotations

import ast
from pathlib import Path
from typing import Any

import pandas as pd

from .base import Step, StepContext, StepError
from ..manifests import structure_id_from_name, write_csv


def _collapse_ranges(indices: list[int], chain_id: str) -> str:
    if not indices:
        return ""
    indices = sorted(set(int(i) for i in indices))
    ranges: list[tuple[int, int]] = []
    start = prev = indices[0]
    for i in indices[1:]:
        if i == prev + 1:
            prev = i
            continue
        ranges.append((start, prev))
        start = prev = i
    ranges.append((start, prev))
    parts = []
    for a, b in ranges:
        if a == b:
            parts.append(f"{chain_id}{a}")
        else:
            parts.append(f"{chain_id}{a}-{chain_id}{b}")
    return ",".join(parts)


class InterfaceEnrichStep(Step):
    name = "interface_enrich"
    stage = "rosetta"
    supports_indices = False

    def expected_total(self, ctx: StepContext) -> int:
        return 1

    def scan_done(self, ctx: StepContext) -> set[int]:
        out_dir = self.output_dir(ctx)
        if (out_dir / "fixed_positions.csv").exists():
            return {0}
        return set()

    def _default_residue_energy_path(self, ctx: StepContext) -> Path:
        run_dir = ctx.out_dir / "output"
        return run_dir / "rosetta_interface" / "residue_energy.csv"

    def run_full(self, ctx: StepContext) -> None:
        residue_energy_path = self.cfg.get("residue_energy_csv")
        if residue_energy_path:
            p = Path(residue_energy_path)
            if not p.is_absolute():
                p = ctx.out_dir / p
        else:
            p = self._default_residue_energy_path(ctx)
        if not p.exists():
            raise StepError(f"residue_energy.csv not found at {p}")

        df = pd.read_csv(p)
        if "binder_energy" not in df.columns:
            raise StepError("residue_energy.csv missing binder_energy column")

        protocol = ctx.input_data.get("protocol")
        binder_chain = str(ctx.input_data.get("binder_chain") or "A")
        if protocol in {"antibody", "vhh"}:
            framework = ctx.input_data.get("framework") or {}
            binder_chain = str(framework.get("heavy_chain") or binder_chain)

        cutoff = float((ctx.input_data.get("filters") or {}).get("rosetta", {}).get("interface_energy_min") or -5.0)

        rows: list[dict[str, Any]] = []
        for _, row in df.iterrows():
            name = str(row.get("pdbname") or Path(str(row.get("pdbpath", ""))).stem)
            structure_id = structure_id_from_name(name)
            try:
                energy_dict = ast.literal_eval(row.get("binder_energy") or "{}")
            except Exception:
                energy_dict = {}
            residues = [int(k) for k, v in energy_dict.items() if float(v) < cutoff]
            rows.append({
                "structure_id": structure_id,
                "pdb_name": name,
                "pdb_path": row.get("pdbpath"),
                "binder_chain": binder_chain,
                "fixed_positions": ",".join(f"{binder_chain}{r}" for r in sorted(set(residues))),
                "fixed_positions_indices": ",".join(str(r) for r in sorted(set(residues))),
                "motif_contig": _collapse_ranges(residues, binder_chain),
                "num_fixed_positions": len(set(residues)),
            })

        if not rows:
            raise StepError("No interface residues found; cannot build fixed positions")

        # Merge per structure_id by taking union across rows
        merged: dict[str, dict[str, Any]] = {}
        for r in rows:
            sid = r["structure_id"]
            cur = merged.setdefault(sid, {**r, "fixed_positions_set": set()})
            cur["fixed_positions_set"].update([x for x in str(r["fixed_positions_indices"]).split(",") if x])
            if not cur.get("pdb_path") and r.get("pdb_path"):
                cur["pdb_path"] = r.get("pdb_path")
        out_rows: list[dict[str, Any]] = []
        for sid, r in merged.items():
            indices = sorted({int(x) for x in r["fixed_positions_set"]})
            out_rows.append({
                "structure_id": sid,
                "pdb_name": r.get("pdb_name"),
                "pdb_path": r.get("pdb_path"),
                "binder_chain": r.get("binder_chain"),
                "fixed_positions": ",".join(f"{r.get('binder_chain')}{i}" for i in indices),
                "fixed_positions_indices": ",".join(str(i) for i in indices),
                "motif_contig": _collapse_ranges(indices, r.get("binder_chain") or "A"),
                "num_fixed_positions": len(indices),
            })

        output_dir = self.output_dir(ctx)
        out_csv = output_dir / "fixed_positions.csv"
        write_csv(out_csv, out_rows, [
            "structure_id",
            "pdb_name",
            "pdb_path",
            "binder_chain",
            "fixed_positions",
            "fixed_positions_indices",
            "motif_contig",
            "num_fixed_positions",
        ])

    def write_manifest(self, ctx: StepContext) -> None:
        # Manifest is the fixed_positions.csv in output_dir
        return
