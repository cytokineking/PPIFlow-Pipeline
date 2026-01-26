from __future__ import annotations

import csv
import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any

import pandas as pd

from .base import Step, StepContext, StepError
from ..logging_utils import log_command_progress, run_command
from ..manifests import build_name_map, extract_design_id, find_metrics_file, structure_id_from_name, write_csv


class ExternalCommandStep(Step):
    supports_indices = False

    def expected_total(self, ctx: StepContext) -> int:
        # Default to number of samples
        sampling = ctx.input_data.get("sampling") or {}
        return int(sampling.get("samples_per_target", 0) or 0)

    def _manifest_has_rows(self, ctx: StepContext) -> bool:
        manifest = self.cfg.get("manifest")
        if not manifest:
            return False
        p = Path(manifest)
        if not p.is_absolute():
            p = ctx.out_dir / p
        if not p.exists():
            return False
        try:
            with p.open("r", newline="") as handle:
                reader = csv.reader(handle)
                rows = list(reader)
            return len(rows) > 1
        except Exception:
            return False

    def scan_done(self, ctx: StepContext) -> set[int]:
        if self._manifest_has_rows(ctx):
            return set(range(self.expected_total(ctx)))
        return set()

    def run_full(self, ctx: StepContext) -> None:
        cmd = self.cfg.get("command")
        if not cmd:
            raise StepError(
                f"Step {self.name} missing command. Edit {self.cfg.get('config_path')} to provide a command."
            )
        if isinstance(cmd, str):
            cmd = [cmd]
        env = os.environ.copy()
        step_label = str(self.cfg.get("name") or self.name)
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
                step_label,
                1,
                1,
                item="command",
                status=status,
                elapsed=time.time() - start,
                log_file=self.cfg.get("_log_file"),
            )


class SeqDesignStep(ExternalCommandStep):
    name = "seq"
    stage = "seq"

    def _default_chain_list(self, ctx: StepContext) -> str:
        protocol = ctx.input_data.get("protocol")
        if protocol == "binder":
            return str(ctx.input_data.get("binder_chain") or "A")
        framework = ctx.input_data.get("framework") or {}
        return str(framework.get("heavy_chain") or "A")

    def _default_mpnn_run(self, ctx: StepContext) -> Path | None:
        tools = ctx.input_data.get("tools") or {}
        if ctx.input_data.get("protocol") == "binder":
            return Path(tools.get("mpnn_run") or tools.get("mpnn_repo", "")).joinpath("protein_mpnn_run.py")
        return Path(tools.get("abmpnn_run") or tools.get("abmpnn_repo", "")).joinpath("protein_mpnn_run.py")

    def _resolve_weights(self, ckpt: str | None) -> tuple[Path | None, str | None]:
        if not ckpt:
            return None, None
        path = Path(ckpt)
        if path.is_file():
            return path.parent, path.stem
        return path, None

    def _write_bias_jsonl(self, out_dir: Path, residues: list[str], weight: float) -> Path:
        out_dir.mkdir(parents=True, exist_ok=True)
        bias = {r: float(weight) for r in residues}
        path = out_dir / "bias_residues.jsonl"
        path.write_text(json.dumps(bias) + "\n")
        return path

    def _build_fixed_positions_map(self, fixed_positions_csv: str | Path) -> dict[str, dict[str, object]]:
        mapping: dict[str, dict[str, object]] = {}
        path = Path(fixed_positions_csv)
        if not path.exists():
            return mapping
        with path.open("r", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                sid = str(row.get("structure_id") or row.get("pdb_name") or "").strip()
                if not sid:
                    continue
                chain = str(row.get("binder_chain") or row.get("chain") or "").strip()
                indices_raw = row.get("fixed_positions_indices") or row.get("fixed_positions") or ""
                indices = [int(x) for x in re.findall(r"\d+", str(indices_raw))]
                entry = {"chain": chain, "indices": indices}
                mapping[sid] = entry
                pname = str(row.get("pdb_name") or "").strip()
                if pname:
                    mapping.setdefault(pname, entry)
        return mapping

    def _resolve_fixed_positions(
        self,
        pdb_path: Path,
        mapping: dict[str, dict[str, object]],
    ) -> dict[str, object] | None:
        stem = pdb_path.stem
        if "__" in stem:
            sid = stem.split("__", 1)[0]
            if sid in mapping:
                return mapping[sid]
        if stem in mapping:
            return mapping[stem]
        parent = pdb_path.parent.name
        if parent in mapping:
            return mapping[parent]
        return None

    def _write_fixed_positions_jsonl(
        self,
        out_dir: Path,
        pdbs: list[Path],
        fixed_positions_csv: str | Path,
        default_chain: str,
    ) -> Path | None:
        mapping = self._build_fixed_positions_map(fixed_positions_csv)
        if not mapping:
            return None
        out_dir.mkdir(parents=True, exist_ok=True)
        jsonl_path = out_dir / "fixed_positions.jsonl"
        count = 0
        with jsonl_path.open("w") as handle:
            for pdb_path in pdbs:
                entry = self._resolve_fixed_positions(pdb_path, mapping)
                if not entry:
                    continue
                indices = [int(x) for x in entry.get("indices") or []]
                if not indices:
                    continue
                chain = str(entry.get("chain") or default_chain or "A")
                payload = {
                    "name": pdb_path.stem,
                    "fixed_positions": {chain: indices},
                }
                handle.write(json.dumps(payload) + "\n")
                count += 1
        if count == 0:
            return None
        return jsonl_path

    def _merge_fastas(self, vanilla_dir: Path, bias_dir: Path, out_dir: Path) -> None:
        out_dir.mkdir(parents=True, exist_ok=True)
        seq_dir = out_dir / "seqs"
        seq_dir.mkdir(parents=True, exist_ok=True)
        vanilla_seq = vanilla_dir / "seqs"
        bias_seq = bias_dir / "seqs"

        vanilla_files = sorted(vanilla_seq.glob("*.fa*")) if vanilla_seq.exists() else []
        bias_files = sorted(bias_seq.glob("*.fa*")) if bias_seq.exists() else []

        if not vanilla_files and not bias_files:
            return

        for fasta in vanilla_files:
            name = fasta.name
            lines = []
            for line in fasta.read_text().splitlines():
                if line.startswith(">"):
                    line = f"{line}|vanilla"
                lines.append(line)
            bias_path = bias_seq / name
            if bias_path.exists():
                for line in bias_path.read_text().splitlines():
                    if line.startswith(">"):
                        line = f"{line}|bias"
                    lines.append(line)
            (seq_dir / name).write_text("\n".join(lines) + "\n")
        # Handle bias-only files
        for fasta in bias_files:
            name = fasta.name
            if (seq_dir / name).exists():
                continue
            lines = []
            for line in fasta.read_text().splitlines():
                if line.startswith(">"):
                    line = f"{line}|bias"
                lines.append(line)
            (seq_dir / name).write_text("\n".join(lines) + "\n")

    def run_full(self, ctx: StepContext) -> None:
        if self.cfg.get("command"):
            return super().run_full(ctx)

        tools = ctx.input_data.get("tools") or {}
        input_dir = self.cfg.get("input_dir")
        if not input_dir:
            raise StepError(f"Step {self.name} missing input_dir")
        input_dir = str((ctx.out_dir / input_dir) if not Path(input_dir).is_absolute() else input_dir)
        if not Path(input_dir).exists():
            raise StepError(f"Input dir not found: {input_dir}")
        out_dir = str(self.output_dir(ctx))

        mpnn_run = self._default_mpnn_run(ctx)
        if not mpnn_run or not mpnn_run.exists():
            raise StepError(
                f"MPNN runner not found. Set tools.mpnn_run/abmpnn_run or provide command in {self.cfg.get('config_path')}"
            )
        ckpt = tools.get("mpnn_ckpt") if ctx.input_data.get("protocol") == "binder" else tools.get("abmpnn_ckpt")

        seq_cfg = (ctx.input_data.get("sequence_design") or {}).get("round1" if self.cfg.get("name") == "seq1" else "round2") or {}
        num_seq = int(seq_cfg.get("num_seq_per_backbone") or 0)
        sampling_temp = float(seq_cfg.get("sampling_temp") or 0.1)
        if num_seq <= 0:
            raise StepError("sequence_design.num_seq_per_backbone must be > 0")
        chain_list = self._default_chain_list(ctx)
        position_list = self.cfg.get("fixed_positions_csv")
        protocol = ctx.input_data.get("protocol")
        use_soluble = bool(seq_cfg.get("use_soluble_ckpt"))
        if protocol == "binder" and use_soluble and self.cfg.get("name") == "seq2":
            ckpt = tools.get("mpnn_ckpt_soluble") or ckpt

        weight_dir, model_name = self._resolve_weights(ckpt)
        chain_arg = " ".join(str(chain_list).replace(",", " ").split())
        seed = str(int((ctx.state.get("runs") or [{}])[0].get("run_seed", 0) or 0))

        def build_base_cmd(out_folder: Path, num_seq_target: int, temp: float, use_soluble_model: bool) -> list[str]:
            cmd = [
                "python",
                str(mpnn_run),
                "--out_folder",
                str(out_folder),
                "--num_seq_per_target",
                str(num_seq_target),
                "--sampling_temp",
                str(temp),
                "--batch_size",
                str(num_seq_target),
                "--seed",
                seed,
            ]
            if weight_dir is not None:
                cmd.extend(["--path_to_model_weights", str(weight_dir)])
            elif use_soluble_model:
                cmd.append("--use_soluble_model")
            else:
                raise StepError("Missing tools.mpnn_ckpt/abmpnn_ckpt for sequence design")
            if model_name:
                cmd.extend(["--model_name", model_name])
            elif protocol != "binder":
                cmd.extend(["--model_name", "abmpnn"])
            return cmd


        step_label = str(self.cfg.get("name") or self.name)

        pdbs = sorted(Path(input_dir).glob("*.pdb"))
        if not pdbs:
            pdbs = sorted(Path(input_dir).rglob("*.pdb"))
        if not pdbs:
            raise StepError(f"No PDBs found for sequence design in {input_dir}")

        # If PDB basenames collide (e.g., partial flow sample0.pdb in multiple subdirs),
        # stage unique names and use those for MPNN to avoid seq output collisions.
        staged_pdbs = pdbs
        stems = [p.stem for p in pdbs]
        should_stage = self.cfg.get("name") == "seq2" or len(set(stems)) != len(stems)
        if should_stage:
            stage_dir = Path(out_dir) / "pdbs"
            stage_dir.mkdir(parents=True, exist_ok=True)
            staged_pdbs = []
            for pdb_path in pdbs:
                parent = pdb_path.parent.name
                unique_stem = f"{parent}__{pdb_path.stem}" if should_stage else pdb_path.stem
                dst = stage_dir / f"{unique_stem}.pdb"
                if not dst.exists():
                    try:
                        os.link(pdb_path, dst)
                    except OSError:
                        shutil.copy2(pdb_path, dst)
                staged_pdbs.append(dst)
            pdbs = staged_pdbs

        fixed_positions_jsonl = None
        if position_list:
            fixed_positions_jsonl = self._write_fixed_positions_jsonl(
                Path(out_dir),
                pdbs,
                position_list,
                str(chain_list).split(",")[0].strip() or "A",
            )
            if not fixed_positions_jsonl:
                print("Warning: fixed_positions_csv provided but could not generate fixed_positions_jsonl; skipping anchors.")

        env = os.environ.copy()

        bias_enabled = bool(protocol == "binder" and self.cfg.get("name") == "seq1" and seq_cfg.get("bias_large_residues"))
        bias_num = int(seq_cfg.get("bias_num") or 0)
        if bias_enabled and bias_num > 0:
            vanilla_num = max(int(num_seq) - bias_num, 0)
            bias_residues = seq_cfg.get("bias_residues") or ["F", "M", "W"]
            if isinstance(bias_residues, str):
                bias_residues = [r.strip() for r in bias_residues.split(",") if r.strip()]
            bias_weight = float(seq_cfg.get("bias_weight") or 0.7)
            tmp_base = self.output_dir(ctx) / "_tmp_seq1"
            vanilla_dir = tmp_base / "vanilla"
            bias_dir = tmp_base / "bias"

            if vanilla_num > 0:
                base_cmd = build_base_cmd(vanilla_dir, vanilla_num, sampling_temp, use_soluble)
                total = len(pdbs)
                for idx, pdb_path in enumerate(pdbs, start=1):
                    cmd = base_cmd + ["--pdb_path", str(pdb_path)]
                    if chain_arg:
                        cmd.extend(["--pdb_path_chains", chain_arg])
                    if fixed_positions_jsonl:
                        cmd.extend(["--fixed_positions_jsonl", str(fixed_positions_jsonl)])
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
                            step_label,
                            idx,
                            total,
                            item=pdb_path.stem,
                            phase="vanilla",
                            status=status,
                            elapsed=time.time() - start,
                            log_file=self.cfg.get("_log_file"),
                        )

            bias_jsonl = self._write_bias_jsonl(bias_dir, bias_residues, bias_weight)
            base_cmd = build_base_cmd(bias_dir, bias_num, sampling_temp, use_soluble)
            base_cmd.extend(["--bias_AA_jsonl", str(bias_jsonl)])
            total = len(pdbs)
            for idx, pdb_path in enumerate(pdbs, start=1):
                cmd = base_cmd + ["--pdb_path", str(pdb_path)]
                if chain_arg:
                    cmd.extend(["--pdb_path_chains", chain_arg])
                if fixed_positions_jsonl:
                    cmd.extend(["--fixed_positions_jsonl", str(fixed_positions_jsonl)])
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
                        step_label,
                        idx,
                        total,
                        item=pdb_path.stem,
                        phase="bias",
                        status=status,
                        elapsed=time.time() - start,
                        log_file=self.cfg.get("_log_file"),
                    )

            self._merge_fastas(vanilla_dir, bias_dir, Path(out_dir))
            shutil.rmtree(tmp_base, ignore_errors=True)
        else:
            base_cmd = build_base_cmd(Path(out_dir), num_seq, sampling_temp, use_soluble)
            total = len(pdbs)
            for idx, pdb_path in enumerate(pdbs, start=1):
                cmd = base_cmd + ["--pdb_path", str(pdb_path)]
                if chain_arg:
                    cmd.extend(["--pdb_path_chains", chain_arg])
                if fixed_positions_jsonl:
                    cmd.extend(["--fixed_positions_jsonl", str(fixed_positions_jsonl)])
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
                        step_label,
                        idx,
                        total,
                        item=pdb_path.stem,
                        status=status,
                        elapsed=time.time() - start,
                        log_file=self.cfg.get("_log_file"),
                    )

    def write_manifest(self, ctx: StepContext) -> None:
        out_dir = self.output_dir(ctx)
        rows = []
        for fp in sorted(out_dir.rglob("*.fa")) + sorted(out_dir.rglob("*.fasta")):
            rows.append({
                "design_id": extract_design_id(fp.stem),
                "structure_id": structure_id_from_name(fp.stem),
                "fasta_path": str(fp),
            })
        if not rows:
            return
        write_csv(self.manifest_path(ctx), rows, ["design_id", "structure_id", "fasta_path"])

    def scan_done(self, ctx: StepContext) -> set[int]:
        done = super().scan_done(ctx)
        if done:
            return done
        out_dir = self.output_dir(ctx)
        seq_dir = out_dir / "seqs"
        fasta_files = list(seq_dir.glob("*.fa*")) if seq_dir.exists() else []
        if not fasta_files:
            fasta_files = list(out_dir.glob("*.fa*"))
        if fasta_files:
            return set(range(self.expected_total(ctx)))
        return set()


class FlowPackerStep(ExternalCommandStep):
    name = "flowpacker"
    stage = "score"

    def write_manifest(self, ctx: StepContext) -> None:
        out_dir = self.output_dir(ctx)
        rows = []
        for fp in sorted(out_dir.rglob("*.pdb")):
            rows.append({
                "design_id": extract_design_id(fp.stem),
                "structure_id": structure_id_from_name(fp.stem),
                "pdb_path": str(fp),
            })
        if not rows:
            return
        write_csv(self.manifest_path(ctx), rows, ["design_id", "structure_id", "pdb_path"])

    def scan_done(self, ctx: StepContext) -> set[int]:
        done = super().scan_done(ctx)
        if done:
            return done
        out_dir = self.output_dir(ctx)
        for sub in [out_dir / "flowpacker_outputs", out_dir / "after_pdbs", out_dir]:
            if sub.exists() and list(sub.rglob("*.pdb")):
                return set(range(self.expected_total(ctx)))
        return set()


class AF3ScoreStep(ExternalCommandStep):
    name = "af3score"
    stage = "score"

    def write_manifest(self, ctx: StepContext) -> None:
        out_dir = self.output_dir(ctx)
        metrics_path = find_metrics_file(out_dir)
        if not metrics_path:
            return
        df = pd.read_csv(metrics_path)
        name_map = build_name_map(out_dir)

        def _get(row, *keys):
            for k in keys:
                if k in row and pd.notna(row[k]):
                    return row[k]
            return None

        def _get_ci(row, key: str):
            lower_map = {str(c).lower(): c for c in row.index}
            col = lower_map.get(key.lower())
            if col is not None and pd.notna(row[col]):
                return row[col]
            return None

        def _get_pairwise_iptm(row, chain_a: str, chain_b: str):
            val = _get_ci(row, f"iptm_{chain_a}_{chain_b}")
            if val is not None:
                return val
            return _get_ci(row, f"iptm_{chain_b}_{chain_a}")

        rows: list[dict[str, Any]] = []
        iptm_global_col: list[float | None] = []
        iptm_binder_target_col: list[float | None] = []
        protocol = ctx.input_data.get("protocol")
        is_antibody = protocol == "antibody"
        for _, row in df.iterrows():
            desc = _get(row, "description", "name", "model", "pdb_name")
            desc = str(desc) if desc is not None else ""
            design_id = extract_design_id(desc)
            iptm_global = _get(row, "iptm", "ipTM", "AF3Score_interchain_iptm", "AF3Score_chain_iptm")
            iptm_binder_target = None
            if is_antibody:
                iptm_ab = _get_pairwise_iptm(row, "A", "B")
                iptm_bc = _get_pairwise_iptm(row, "B", "C")
                if iptm_ab is None or iptm_bc is None:
                    raise StepError(
                        "AF3Score metrics missing pairwise iptm columns required for antibody scoring "
                        "(need iptm_A_B and iptm_B_C or reversed ordering)."
                    )
                iptm_binder_target = (float(iptm_ab) + float(iptm_bc)) / 2.0
                iptm = iptm_binder_target
            else:
                iptm = iptm_global
            ptm = _get(row, "ptm", "pTM", "ptm_A", "ptm_B")
            pdb_path = name_map.get(desc)
            iptm_global_col.append(float(iptm_global) if iptm_global is not None else None)
            iptm_binder_target_col.append(
                float(iptm_binder_target) if iptm_binder_target is not None else None
            )
            rows.append({
                "design_id": design_id,
                "structure_id": structure_id_from_name(desc),
                "iptm": float(iptm) if iptm is not None else None,
                "iptm_binder_target": (
                    float(iptm_binder_target) if iptm_binder_target is not None else None
                ),
                "iptm_global": float(iptm_global) if iptm_global is not None else None,
                "ptm": float(ptm) if ptm is not None else None,
                "pdb_path": str(pdb_path) if pdb_path else None,
            })

        # Write derived metrics with binder-target ipTM if needed
        df = df.copy()
        df["iptm_global"] = iptm_global_col
        df["iptm_binder_target"] = iptm_binder_target_col
        metrics_ppiflow = out_dir / "metrics_ppiflow.csv"
        try:
            df.to_csv(metrics_ppiflow, index=False)
        except Exception:
            pass

        if not rows:
            return

        filters = (ctx.input_data.get("filters") or {}).get("af3score") or {}
        if self.cfg.get("name") == "af3score2":
            iptm_min = float((filters.get("round2") or {}).get("iptm_min") or 0)
            ptm_min = float((filters.get("round2") or {}).get("ptm_min") or 0)
        else:
            iptm_min = float((filters.get("round1") or {}).get("iptm_min") or 0)
            ptm_val = (filters.get("round1") or {}).get("ptm_min")
            ptm_min = float(ptm_val or 0) if ptm_val is not None else 0.0

        for r in rows:
            iptm_ok = r.get("iptm") is not None and float(r.get("iptm") or 0) >= iptm_min
            ptm_ok = True
            if ptm_min:
                ptm_ok = r.get("ptm") is not None and float(r.get("ptm") or 0) >= ptm_min
            r["passed_filter"] = bool(iptm_ok and ptm_ok)

        top_k = None
        if self.cfg.get("name") == "af3score2":
            top_k = (filters.get("round2") or {}).get("top_k")
        else:
            top_k = (filters.get("round1") or {}).get("top_k")
        try:
            top_k = int(top_k) if top_k is not None else None
        except Exception:
            top_k = None

        filtered_rows = [r for r in rows if r.get("passed_filter")]
        if top_k is not None and top_k > 0:
            filtered_rows = sorted(
                filtered_rows,
                key=lambda r: (r.get("iptm") is not None, float(r.get("iptm") or 0)),
                reverse=True,
            )
            keep_ids = set()
            for r in filtered_rows[:top_k]:
                keep_ids.add((r.get("design_id"), r.get("structure_id")))
            for r in rows:
                r["passed_top_k"] = (r.get("design_id"), r.get("structure_id")) in keep_ids
        else:
            for r in rows:
                r["passed_top_k"] = r.get("passed_filter")

        # Write filtered PDBs for downstream steps
        filtered_dir = out_dir / "filtered_pdbs"
        filtered_dir.mkdir(parents=True, exist_ok=True)
        for r in rows:
            if not r.get("passed_filter"):
                continue
            if top_k is not None and not r.get("passed_top_k"):
                continue
            pdb_path = r.get("pdb_path")
            if not pdb_path:
                continue
            src = Path(str(pdb_path))
            if not src.exists():
                continue
            dst = filtered_dir / src.name
            if dst.exists():
                continue
            try:
                os.link(src, dst)
            except Exception:
                try:
                    dst.write_bytes(src.read_bytes())
                except Exception:
                    pass

        write_csv(
            self.manifest_path(ctx),
            rows,
            [
                "design_id",
                "structure_id",
                "iptm",
                "iptm_binder_target",
                "iptm_global",
                "ptm",
                "pdb_path",
                "passed_filter",
                "passed_top_k",
            ],
        )

    def scan_done(self, ctx: StepContext) -> set[int]:
        done = super().scan_done(ctx)
        if done:
            return done
        out_dir = self.output_dir(ctx)
        metrics = out_dir / "metrics.csv"
        if metrics.exists():
            return set(range(self.expected_total(ctx)))
        return set()


class RosettaInterfaceStep(ExternalCommandStep):
    name = "rosetta_interface"
    stage = "rosetta"

    def write_manifest(self, ctx: StepContext) -> None:
        out_dir = self.output_dir(ctx)
        residue_csv = out_dir / "residue_energy.csv"
        if not residue_csv.exists():
            return
        try:
            df = pd.read_csv(residue_csv)
        except Exception:
            return
        rows = []
        for _, row in df.iterrows():
            name = str(row.get("pdbname") or Path(str(row.get("pdbpath", ""))).stem)
            rows.append({
                "design_id": extract_design_id(name),
                "structure_id": structure_id_from_name(name),
                "residue_energy_csv": str(residue_csv),
            })
        if not rows:
            return
        write_csv(self.manifest_path(ctx), rows, ["design_id", "structure_id", "residue_energy_csv"])

    def scan_done(self, ctx: StepContext) -> set[int]:
        done = super().scan_done(ctx)
        if done:
            return done
        out_dir = self.output_dir(ctx)
        if (out_dir / "residue_energy.csv").exists():
            return set(range(self.expected_total(ctx)))
        return set()


class RelaxStep(ExternalCommandStep):
    name = "relax"
    stage = "rosetta"

    def write_manifest(self, ctx: StepContext) -> None:
        out_dir = self.output_dir(ctx)
        rows = []
        for fp in sorted(out_dir.rglob("*.pdb")):
            rows.append({
                "design_id": extract_design_id(fp.stem),
                "structure_id": structure_id_from_name(fp.stem),
                "pdb_path": str(fp),
            })
        if not rows:
            return
        write_csv(self.manifest_path(ctx), rows, ["design_id", "structure_id", "pdb_path"])

    def scan_done(self, ctx: StepContext) -> set[int]:
        done = super().scan_done(ctx)
        if done:
            return done
        out_dir = self.output_dir(ctx)
        if list(out_dir.rglob("*.pdb")):
            return set(range(self.expected_total(ctx)))
        return set()


class AF3RefoldStep(ExternalCommandStep):
    name = "af3_refold"
    stage = "score"

    def write_manifest(self, ctx: StepContext) -> None:
        out_dir = self.output_dir(ctx)
        metrics_path = find_metrics_file(out_dir)
        if not metrics_path:
            return
        try:
            df = pd.read_csv(metrics_path)
        except Exception:
            return
        name_map = build_name_map(out_dir / "pdbs")

        def _get(row, *keys):
            for k in keys:
                if k in row and pd.notna(row[k]):
                    return row[k]
            return None

        rows: list[dict[str, Any]] = []
        for _, row in df.iterrows():
            desc = _get(row, "description", "name", "model", "pdb_name")
            desc = str(desc) if desc is not None else ""
            if not desc:
                continue
            iptm = _get(row, "iptm", "ipTM", "AF3Score_interchain_iptm", "AF3Score_chain_iptm")
            ptm = _get(row, "ptm", "pTM", "ptm_A", "ptm_B", "AF3Score_chain_ptm")
            pdb_path = name_map.get(desc)
            rows.append({
                "design_id": extract_design_id(desc),
                "structure_id": structure_id_from_name(desc),
                "iptm": float(iptm) if iptm is not None else None,
                "ptm": float(ptm) if ptm is not None else None,
                "pdb_path": str(pdb_path) if pdb_path else None,
            })

        if rows:
            write_csv(self.manifest_path(ctx), rows, ["design_id", "structure_id", "iptm", "ptm", "pdb_path"])

    def run_full(self, ctx: StepContext) -> None:
        cmd = self.cfg.get("command")
        if not cmd:
            return
        return super().run_full(ctx)

    def scan_done(self, ctx: StepContext) -> set[int]:
        done = super().scan_done(ctx)
        if done:
            return done
        out_dir = self.output_dir(ctx)
        metrics = out_dir / "metrics.csv"
        if metrics.exists():
            return set(range(self.expected_total(ctx)))
        return set()
