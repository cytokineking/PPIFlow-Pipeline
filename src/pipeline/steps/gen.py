from __future__ import annotations

import os
import subprocess
import time
from pathlib import Path
from typing import Any

from .base import Step, StepContext, StepError
from ..logging_utils import log_command_progress, run_command
from ..io import repo_root, write_json
from ..manifests import extract_design_id, structure_id_from_name, write_csv


def _resolve_repo_path(path: str | Path) -> Path:
    p = Path(path)
    if p.is_absolute():
        return p
    return repo_root() / p


class GenStep(Step):
    name = "gen"
    stage = "gen"
    supports_indices = True

    def expected_total(self, ctx: StepContext) -> int:
        sampling = ctx.input_data.get("sampling") or {}
        n = int(sampling.get("samples_per_target", 0) or 0)
        if n <= 0:
            raise StepError("sampling.samples_per_target must be > 0")
        return n

    def _parse_existing(self, out_dir: Path, name: str) -> set[int]:
        done: set[int] = set()
        if not out_dir.exists():
            return done
        for fp in out_dir.glob(f"{name}_*.pdb"):
            try:
                stem = fp.stem
                if not stem.startswith(f"{name}_"):
                    continue
                suffix = stem[len(name) + 1 :]
                if suffix.isdigit():
                    done.add(int(suffix))
            except Exception:
                continue
        return done

    def scan_done(self, ctx: StepContext) -> set[int]:
        out_dir = self.output_dir(ctx)
        name = ctx.input_data.get("name", "")
        done = self._parse_existing(out_dir, name)
        expected = self.expected_total(ctx)
        if done and len(done) == expected:
            return done
        # Fallback: include non-numeric stems when count matches expected.
        if expected > 0:
            pdbs = sorted(out_dir.glob(f"{name}_*.pdb"))
            if len(pdbs) == expected:
                return set(range(expected))
        return done

    def _length_schedule_path(self, ctx: StepContext, out_dir: Path) -> Path:
        path = out_dir / "length_schedule.json"
        return path

    def _ensure_length_schedule(self, ctx: StepContext, out_dir: Path) -> Path | None:
        protocol = ctx.input_data.get("protocol")
        if protocol != "binder":
            return None
        binder = ctx.input_data.get("binder") or {}
        length_spec = binder.get("length")
        if not length_spec:
            return None
        if "-" in str(length_spec):
            parts = str(length_spec).split("-")
            min_len, max_len = int(parts[0]), int(parts[1])
        else:
            min_len = max_len = int(length_spec)
        schedule_path = self._length_schedule_path(ctx, out_dir)
        if schedule_path.exists():
            return schedule_path
        # Generate deterministic schedule
        sampling = ctx.input_data.get("sampling") or {}
        n = int(sampling.get("samples_per_target", 0) or 0)
        seed = int((ctx.state.get("runs") or [{}])[0].get("run_seed", 0) or 0)
        try:
            import numpy as np

            rng = np.random.default_rng(seed=seed)
            lengths = rng.integers(low=min_len, high=max_len, size=(n,))
            schedule = [int(x) for x in lengths]
        except Exception:
            schedule = [int(min_len)] * n
        write_json(schedule_path, {"lengths": schedule}, indent=2)
        return schedule_path

    def run_indices(self, ctx: StepContext, indices: list[int]) -> None:
        protocol = ctx.input_data.get("protocol")
        out_dir = self.output_dir(ctx)
        name = ctx.input_data.get("name")
        tools = ctx.input_data.get("tools") or {}
        config_path = tools.get("ppiflow_config")

        if protocol == "binder":
            script = _resolve_repo_path("src/entrypoints/sample_binder.py")
            config_default = tools.get("ppiflow_binder_config") or str(
                _resolve_repo_path("src/configs/inference_binder.yaml")
            )
            config_path = config_path or config_default
            config_path = str(_resolve_repo_path(config_path))
            binder = ctx.input_data.get("binder") or {}
            length_spec = binder.get("length")
            if not length_spec:
                raise StepError("binder.length missing")
            if "-" in str(length_spec):
                min_len, max_len = [int(x) for x in str(length_spec).split("-")]
            else:
                min_len = max_len = int(length_spec)

            target = ctx.input_data.get("target") or {}
            target_chains = target.get("chains") or []
            target_chain = target_chains[0] if target_chains else None
            binder_chain = ctx.input_data.get("binder_chain") or "A"
            hotspots = target.get("hotspots")
            hotspots_file = target.get("hotspots_file")
            ckpt = tools.get("ppiflow_ckpt")
            if not ckpt:
                raise StepError("tools.ppiflow_ckpt is required for binder generation")

            cmd = [
                "python",
                str(script),
                "--input_pdb",
                str(target.get("pdb")),
                "--target_chain",
                str(target_chain),
                "--binder_chain",
                str(binder_chain),
                "--config",
                str(config_path),
                "--samples_min_length",
                str(min_len),
                "--samples_max_length",
                str(max_len),
                "--samples_per_target",
                str(self.expected_total(ctx)),
                "--model_weights",
                str(ckpt),
                "--output_dir",
                str(out_dir),
                "--name",
                str(name),
                "--seed",
                str(int((ctx.state.get("runs") or [{}])[0].get("run_seed", 0) or 0)),
                "--sample_ids",
                ",".join(str(i) for i in indices),
            ]
            if hotspots_file:
                cmd.extend(["--hotspots_file", str(hotspots_file)])
            elif hotspots:
                if isinstance(hotspots, list):
                    cmd.extend(["--specified_hotspots", ",".join(hotspots)])
                else:
                    cmd.extend(["--specified_hotspots", str(hotspots)])
            schedule_path = self._ensure_length_schedule(ctx, out_dir)
            if schedule_path is not None:
                cmd.extend(["--length_schedule_path", str(schedule_path)])
        else:
            script = _resolve_repo_path("src/entrypoints/sample_antibody_nanobody.py")
            config_default = tools.get("ppiflow_antibody_config") or str(
                _resolve_repo_path("src/configs/inference_nanobody.yaml")
            )
            config_path = config_path or config_default
            config_path = str(_resolve_repo_path(config_path))
            target = ctx.input_data.get("target") or {}
            framework = ctx.input_data.get("framework") or {}
            ckpt = tools.get("ppiflow_ckpt")
            if not ckpt:
                raise StepError("tools.ppiflow_ckpt is required for antibody/vhh generation")
            cmd = [
                "python",
                str(script),
                "--antigen_pdb",
                str(target.get("pdb")),
                "--framework_pdb",
                str(framework.get("pdb")),
                "--antigen_chain",
                ",".join(target.get("chains") or []),
                "--heavy_chain",
                str(framework.get("heavy_chain")),
                "--cdr_length",
                str(framework.get("cdr_length")),
                "--samples_per_target",
                str(self.expected_total(ctx)),
                "--config",
                str(config_path),
                "--model_weights",
                str(ckpt),
                "--output_dir",
                str(out_dir),
                "--name",
                str(name),
                "--seed",
                str(int((ctx.state.get("runs") or [{}])[0].get("run_seed", 0) or 0)),
                "--sample_ids",
                ",".join(str(i) for i in indices),
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

        env = os.environ.copy()
        start = time.time()
        status = "OK"
        try:
            run_command(
                cmd,
                env=env,
                log_file=self.cfg.get("_log_file"),
                verbose=bool(self.cfg.get("_verbose")),
            )
        except Exception:
            status = "FAILED"
            raise
        finally:
            log_command_progress(
                str(self.cfg.get("name") or self.name),
                1,
                1,
                item="command",
                status=status,
                elapsed=time.time() - start,
                log_file=self.cfg.get("_log_file"),
            )

    def write_manifest(self, ctx: StepContext) -> None:
        out_dir = self.output_dir(ctx)
        name = ctx.input_data.get("name", "design")
        rows = []
        for fp in sorted(out_dir.glob(f"{name}_*.pdb")):
            did = extract_design_id(fp.stem)
            rows.append({
                "design_id": did,
                "structure_id": structure_id_from_name(fp.stem),
                "pdb_path": str(fp),
                "parent_id": None,
            })
        if not rows:
            return
        manifest = self.manifest_path(ctx)
        write_csv(manifest, rows, ["design_id", "structure_id", "pdb_path", "parent_id"])
