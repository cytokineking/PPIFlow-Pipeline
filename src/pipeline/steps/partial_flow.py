from __future__ import annotations

import os
import subprocess
import time
from pathlib import Path
from typing import Any

import pandas as pd

from .base import Step, StepContext, StepError
from ..logging_utils import log_command_progress, run_command
from ..io import repo_root
from ..manifests import extract_design_id, structure_id_from_name, write_csv


def _resolve_repo_path(path: str | Path) -> Path:
    p = Path(path)
    if p.is_absolute():
        return p
    return repo_root() / p


class PartialFlowStep(Step):
    name = "partial"
    stage = "partial"
    supports_indices = True

    def expected_total(self, ctx: StepContext) -> int:
        rows = self._load_fixed_positions(ctx)
        if not rows:
            raise StepError("No fixed positions found for partial flow")
        return len(rows)

    def _load_fixed_positions(self, ctx: StepContext) -> list[dict[str, Any]]:
        fixed_positions_csv = self.cfg.get("fixed_positions_csv")
        if not fixed_positions_csv:
            return []
        p = Path(fixed_positions_csv)
        if not p.is_absolute():
            p = ctx.out_dir / p
        if not p.exists():
            return []
        df = pd.read_csv(p)
        rows = df.to_dict(orient="records")
        return rows

    def scan_done(self, ctx: StepContext) -> set[int]:
        out_dir = self.output_dir(ctx)
        rows = self._load_fixed_positions(ctx)
        done: set[int] = set()
        for idx, row in enumerate(rows):
            sid = str(row.get("structure_id") or row.get("pdb_name") or idx)
            sub = out_dir / sid
            if sub.exists():
                if list(sub.rglob("sample*.pdb")) or list(sub.rglob("sample*_*.pdb")):
                    done.add(idx)
        return done

    def run_indices(self, ctx: StepContext, indices: list[int]) -> None:
        rows = self._load_fixed_positions(ctx)
        if not rows:
            raise StepError("fixed_positions_csv missing or empty")
        protocol = ctx.input_data.get("protocol")
        out_dir = self.output_dir(ctx)
        name = ctx.input_data.get("name")
        tools = ctx.input_data.get("tools") or {}
        partial_cfg = ctx.input_data.get("partial") or {}
        samples_per_target = int(partial_cfg.get("samples_per_target") or 8)

        if protocol == "binder":
            script = _resolve_repo_path("src/entrypoints/sample_binder_partial.py")
            config_path = tools.get("ppiflow_binder_partial_config") or str(
                _resolve_repo_path("src/configs/inference_binder_partial.yaml")
            )
            config_path = str(_resolve_repo_path(config_path))
            target = ctx.input_data.get("target") or {}
            binder_chain = ctx.input_data.get("binder_chain") or "A"
            target_chain = target.get("chains", [None])[0]
            ckpt = tools.get("ppiflow_ckpt")
            if not ckpt:
                raise StepError("tools.ppiflow_ckpt is required for binder partial flow")
            total = len(indices)
            for pos, idx in enumerate(indices, start=1):
                row = rows[idx]
                sid = str(row.get("structure_id") or row.get("pdb_name") or idx)
                pdb_path = row.get("pdb_path") or target.get("pdb")
                motif_contig = row.get("motif_contig")
                out_sub = out_dir / sid
                cmd = [
                    "python",
                    str(script),
                    "--input_pdb",
                    str(pdb_path),
                    "--target_chain",
                    str(target_chain),
                    "--binder_chain",
                    str(binder_chain),
                    "--config",
                    str(config_path),
                    "--samples_per_target",
                    str(samples_per_target),
                    "--model_weights",
                    str(ckpt),
                    "--output_dir",
                    str(out_sub),
                    "--name",
                    str(name),
                    "--seed",
                    str(int((ctx.state.get("runs") or [{}])[0].get("run_seed", 0) or 0)),
                ]
                if motif_contig:
                    cmd.extend(["--motif_contig", str(motif_contig)])
                hotspots = target.get("hotspots")
                hotspots_file = target.get("hotspots_file")
                if hotspots_file:
                    cmd.extend(["--hotspots_file", str(hotspots_file)])
                elif hotspots:
                    if isinstance(hotspots, list):
                        cmd.extend(["--specified_hotspots", ",".join(hotspots)])
                    else:
                        cmd.extend(["--specified_hotspots", str(hotspots)])
                start = time.time()
                status = "OK"
                try:
                    run_command(
                        cmd,
                        env=os.environ.copy(),
                        log_file=self.cfg.get("_log_file"),
                        verbose=bool(self.cfg.get("_verbose")),
                    )
                except Exception:
                    status = "FAILED"
                    raise
                finally:
                    log_command_progress(
                        str(self.cfg.get("name") or self.name),
                        pos,
                        total,
                        item=sid,
                        status=status,
                        elapsed=time.time() - start,
                        log_file=self.cfg.get("_log_file"),
                    )
        else:
            script = _resolve_repo_path("src/entrypoints/sample_antibody_nanobody_partial.py")
            config_path = tools.get("ppiflow_antibody_partial_config") or str(
                _resolve_repo_path("src/configs/inference_nanobody.yaml")
            )
            config_path = str(_resolve_repo_path(config_path))
            target = ctx.input_data.get("target") or {}
            framework = ctx.input_data.get("framework") or {}
            ckpt = tools.get("ppiflow_ckpt")
            if not ckpt:
                raise StepError("tools.ppiflow_ckpt is required for antibody/vhh partial flow")
            total = len(indices)
            for pos, idx in enumerate(indices, start=1):
                row = rows[idx]
                sid = str(row.get("structure_id") or row.get("pdb_name") or idx)
                pdb_path = row.get("pdb_path")
                if not pdb_path:
                    raise StepError("fixed_positions.csv missing pdb_path for partial flow")
                fixed_positions = row.get("fixed_positions")
                if not fixed_positions:
                    raise StepError("fixed_positions.csv missing fixed_positions")
                out_sub = out_dir / sid
                cmd = [
                    "python",
                    str(script),
                    "--complex_pdb",
                    str(pdb_path),
                    "--fixed_positions",
                    str(fixed_positions),
                    "--antigen_chain",
                    ",".join(target.get("chains") or []),
                    "--heavy_chain",
                    str(framework.get("heavy_chain")),
                    "--samples_per_target",
                    str(samples_per_target),
                    "--config",
                    str(config_path),
                    "--model_weights",
                    str(ckpt),
                    "--output_dir",
                    str(out_sub),
                    "--name",
                    str(name),
                    "--start_t",
                    str(ctx.input_data.get("partial", {}).get("start_t", 0.6)),
                    "--seed",
                    str(int((ctx.state.get("runs") or [{}])[0].get("run_seed", 0) or 0)),
                ]
                light_chain = framework.get("light_chain")
                if light_chain:
                    cmd.extend(["--light_chain", str(light_chain)])
                hotspots = target.get("hotspots")
                hotspots_file = target.get("hotspots_file")
                if hotspots_file:
                    cmd.extend(["--hotspots_file", str(hotspots_file)])
                elif hotspots:
                    if isinstance(hotspots, list):
                        cmd.extend(["--specified_hotspots", ",".join(hotspots)])
                    else:
                        cmd.extend(["--specified_hotspots", str(hotspots)])
                start = time.time()
                status = "OK"
                try:
                    run_command(
                        cmd,
                        env=os.environ.copy(),
                        log_file=self.cfg.get("_log_file"),
                        verbose=bool(self.cfg.get("_verbose")),
                    )
                except Exception:
                    status = "FAILED"
                    raise
                finally:
                    log_command_progress(
                        str(self.cfg.get("name") or self.name),
                        pos,
                        total,
                        item=sid,
                        status=status,
                        elapsed=time.time() - start,
                        log_file=self.cfg.get("_log_file"),
                    )

    def write_manifest(self, ctx: StepContext) -> None:
        out_dir = self.output_dir(ctx)
        rows = []
        for fp in out_dir.rglob("sample*.pdb"):
            stem = fp.stem
            did = extract_design_id(stem) or extract_design_id(fp.parent.name)
            rows.append({
                "design_id": did,
                "structure_id": structure_id_from_name(fp.parent.name),
                "pdb_path": str(fp),
            })
        if not rows:
            return
        write_csv(self.manifest_path(ctx), rows, ["design_id", "structure_id", "pdb_path"])
