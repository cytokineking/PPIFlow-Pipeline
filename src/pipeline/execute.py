from __future__ import annotations

import os
import sys
import time
import traceback
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
import csv
from typing import Iterable

from .heartbeat import HeartbeatReporter
from .io import read_yaml
from .state import collect_tool_versions, init_or_update_state, load_state, sha256_json
from .steps import STEP_GROUPS, STEP_ORDER, STEP_REGISTRY
from .steps.base import StepContext, StepError


def _resolve_steps(steps_arg: str | None) -> list[str]:
    if not steps_arg or steps_arg == "all":
        return STEP_ORDER
    tokens = [s.strip() for s in steps_arg.split(",") if s.strip()]
    resolved: list[str] = []
    for tok in tokens:
        if tok in STEP_GROUPS:
            for s in STEP_GROUPS[tok]:
                if s not in resolved:
                    resolved.append(s)
        elif tok in STEP_ORDER:
            if tok not in resolved:
                resolved.append(tok)
        else:
            raise StepError(f"Unknown step group or step: {tok}")
    return resolved


def _console(msg: str) -> None:
    print(msg, file=sys.__stdout__, flush=True)


def _format_duration(seconds: float) -> str:
    seconds = int(seconds)
    mins, sec = divmod(seconds, 60)
    hrs, mins = divmod(mins, 60)
    if hrs:
        return f"{hrs:02d}h{mins:02d}m{sec:02d}s"
    return f"{mins:02d}m{sec:02d}s"


def _count_csv_rows(path: Path) -> int | None:
    if not path.exists():
        return None
    try:
        with path.open("r", newline="") as handle:
            reader = csv.reader(handle)
            rows = list(reader)
        if not rows:
            return 0
        # assume header
        return max(len(rows) - 1, 0)
    except Exception:
        return None


def _tail_lines(path: Path, max_lines: int = 20) -> list[str]:
    if not path.exists():
        return []
    try:
        with path.open("r") as handle:
            lines = handle.readlines()
        if len(lines) <= max_lines:
            return [line.rstrip("\n") for line in lines]
        return [line.rstrip("\n") for line in lines[-max_lines:]]
    except Exception:
        return []


def _mirror_log_prefixes(log_file: Path, prefixes: list[str]) -> int:
    if not log_file.exists():
        return 0
    count = 0
    try:
        with log_file.open("r") as handle:
            for line in handle:
                for prefix in prefixes:
                    if line.startswith(prefix):
                        _console(line.rstrip("\n"))
                        count += 1
                        break
    except Exception:
        return count
    return count


def _step_params(step_name: str, input_data: dict) -> list[str]:
    params: list[str] = []
    if step_name in {"af3score1", "af3score2"}:
        filters = (input_data.get("filters") or {}).get("af3score") or {}
        if step_name == "af3score1":
            iptm = (filters.get("round1") or {}).get("iptm_min")
            if iptm is not None:
                params.append(f"iptm_min={iptm}")
            ptm = (filters.get("round1") or {}).get("ptm_min")
            if ptm is not None:
                params.append(f"ptm_min={ptm}")
            top_k = (filters.get("round1") or {}).get("top_k")
            if top_k:
                params.append(f"top_k={top_k}")
        else:
            r2 = filters.get("round2") or {}
            if r2.get("iptm_min") is not None:
                params.append(f"iptm_min={r2.get('iptm_min')}")
            if r2.get("ptm_min") is not None:
                params.append(f"ptm_min={r2.get('ptm_min')}")
            if r2.get("top_k"):
                params.append(f"top_k={r2.get('top_k')}")
    if step_name in {"seq1", "seq2"}:
        seq_cfg = (input_data.get("sequence_design") or {}).get("round1" if step_name == "seq1" else "round2") or {}
        if seq_cfg.get("num_seq_per_backbone") is not None:
            params.append(f"num_seq={seq_cfg.get('num_seq_per_backbone')}")
        if seq_cfg.get("sampling_temp") is not None:
            params.append(f"temp={seq_cfg.get('sampling_temp')}")
        if step_name == "seq1" and seq_cfg.get("bias_large_residues"):
            residues = seq_cfg.get("bias_residues") or "F,M,W"
            params.append(f"bias_large_residues={residues}")
    if step_name == "partial":
        start_t = (input_data.get("partial") or {}).get("start_t")
        if start_t is not None:
            params.append(f"start_t={start_t}")
        samples = (input_data.get("partial") or {}).get("samples_per_target")
        if samples is not None:
            params.append(f"samples={samples}")
    if step_name == "rosetta_interface":
        rosetta = (input_data.get("filters") or {}).get("rosetta") or {}
        if rosetta.get("interface_energy_min") is not None:
            params.append(f"interface_energy_min={rosetta.get('interface_energy_min')}")
    if step_name == "af3_refold":
        refold = (input_data.get("filters") or {}).get("af3_refold") or {}
        if refold.get("iptm_min") is not None:
            params.append(f"iptm_min={refold.get('iptm_min')}")
        if refold.get("ptm_min") is not None:
            params.append(f"ptm_min={refold.get('ptm_min')}")
        if refold.get("dockq_min") is not None:
            params.append(f"dockq_min={refold.get('dockq_min')}")
    if step_name == "dockq":
        dockq = (input_data.get("filters") or {}).get("dockq") or {}
        if dockq.get("min") is not None:
            params.append(f"dockq_min={dockq.get('min')}")
    if step_name == "rank":
        ranking = input_data.get("ranking") or {}
        if ranking.get("top_k") is not None:
            params.append(f"top_k={ranking.get('top_k')}")
    return params


def execute_pipeline(args) -> None:
    out_dir = Path(args.output).resolve()
    steps_yaml = out_dir / "steps.yaml"
    if not steps_yaml.exists():
        raise StepError(f"steps.yaml not found at {steps_yaml}")

    steps_data = read_yaml(steps_yaml)
    steps_list = steps_data.get("steps") if isinstance(steps_data, dict) else None
    if not steps_list:
        raise StepError("steps.yaml is empty or invalid")

    input_json = out_dir / "pipeline_input.json"
    if not input_json.exists():
        raise StepError("pipeline_input.json missing; run configure first")

    from .io import read_json

    input_data = read_json(input_json)
    tool_versions = collect_tool_versions(input_data.get("tools") or {})
    sampling = input_data.get("sampling") or {}
    target_n = int(sampling.get("samples_per_target", 0) or 0)
    state = init_or_update_state(
        out_dir=out_dir,
        input_sha256=sha256_json(input_data),
        tool_versions=tool_versions,
        target_n=target_n,
        seeds=(sampling.get("seeds") or None),
    )

    run_id = int((state.get("runs") or [{"run_id": 0}])[0].get("run_id", 0))
    run_seed = int((state.get("runs") or [{"run_seed": 0}])[0].get("run_seed", 0) or 0)

    os.environ["PPIFLOW_PROTOCOL"] = str(input_data.get("protocol"))
    os.environ["PPIFLOW_RUN_ID"] = str(run_id)
    os.environ["PPIFLOW_SEED"] = str(run_seed)
    if getattr(args, "continue_on_error", False):
        options = input_data.setdefault("options", {})
        options.setdefault("continue_on_item_error", True)

    hb = HeartbeatReporter.from_env(out_dir)

    rank = int(os.environ.get("RANK", "0"))
    world_size = int(os.environ.get("WORLD_SIZE", "1"))
    local_rank = int(os.environ.get("LOCAL_RANK", "0"))

    ctx = StepContext(
        out_dir=out_dir,
        input_data=input_data,
        state=state,
        run_id=run_id,
        rank=rank,
        world_size=world_size,
        local_rank=local_rank,
        reuse=bool(args.reuse),
        heartbeat=hb,
    )

    selected_steps = set(_resolve_steps(args.steps))
    ordered_steps: list[dict] = []
    for step_info in steps_list:
        step_name = step_info.get("name")
        if step_name in selected_steps:
            ordered_steps.append(step_info)

    total_steps = len(ordered_steps)
    _console("=" * 78)
    _console("PIPELINE START")
    _console(f"output: {out_dir}")
    _console(f"steps: {total_steps}")
    _console(f"verbose: {bool(getattr(args, 'verbose', False))}")
    _console("=" * 78)

    log_dir = out_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    manifest_counts: dict[str, int] = {}

    for idx, step_info in enumerate(ordered_steps, start=1):
        step_name = step_info.get("name")
        config_rel = step_info.get("config_file")
        if not step_name or not config_rel:
            continue
        cfg_path = out_dir / config_rel
        cfg = read_yaml(cfg_path)
        cfg["config_path"] = str(cfg_path)

        StepCls = STEP_REGISTRY.get(step_name)
        if StepCls is None:
            raise StepError(f"Unknown step: {step_name}")
        step = StepCls(cfg)

        log_file = log_dir / f"{idx:02d}_{step.name}.log"
        cfg["_log_file"] = str(log_file)
        cfg["_verbose"] = bool(getattr(args, "verbose", False))

        os.environ["PPIFLOW_STAGE"] = str(step.stage)
        os.environ["PPIFLOW_STEP"] = str(step.name)
        if hb:
            hb.touch(extra={"stage_transition": step.stage, "step": step.name}, force=True)

        _console("=" * 78)
        _console(f"STEP {idx}/{total_steps} | {step.name} ({step.stage})")
        _console(f"config: {cfg_path}")
        _console(f"log: {log_file}")
        params = _step_params(step_name, input_data)
        if params:
            _console(f"params: {', '.join(params)}")
        _console("-" * 78)

        start_time = time.time()
        failed = False
        err_text = None
        pending_exc = None
        log_handle = None
        try:
            if not getattr(args, "verbose", False):
                log_handle = log_file.open("a")
                with redirect_stdout(log_handle), redirect_stderr(log_handle):
                    step.run(ctx)
            else:
                step.run(ctx)
        except Exception:
            failed = True
            err_text = traceback.format_exc()
            pending_exc = sys.exc_info()
            errors_path = out_dir / "errors.log"
            with errors_path.open("a") as handle:
                handle.write(f"[{step.name}] failed\n")
                handle.write(err_text)
                handle.write("\n")
            if hb:
                hb.touch(extra={"step_error": step.name}, force=True)
        finally:
            if log_handle:
                log_handle.close()

        elapsed = _format_duration(time.time() - start_time)
        input_count = None
        input_dir = cfg.get("input_dir")
        if input_dir:
            input_dir_path = Path(input_dir)
            if not input_dir_path.is_absolute():
                input_dir_path = out_dir / input_dir_path
            input_count = manifest_counts.get(str(input_dir_path.resolve()))

        output_count = None
        output_dir = cfg.get("output_dir")
        if output_dir:
            output_dir_path = Path(output_dir)
            if not output_dir_path.is_absolute():
                output_dir_path = out_dir / output_dir_path
        else:
            output_dir_path = None

        if cfg.get("manifest"):
            manifest_path = Path(cfg["manifest"])
            if not manifest_path.is_absolute():
                manifest_path = out_dir / manifest_path
            output_count = _count_csv_rows(manifest_path)

        if output_dir_path is not None and output_count is not None:
            manifest_counts[str(output_dir_path.resolve())] = output_count

        status = "FAILED" if failed else "OK"
        _console(
            f"RESULT: {status} | elapsed={elapsed} | in={input_count if input_count is not None else 'n/a'} | out={output_count if output_count is not None else 'n/a'}"
        )
        if failed:
            _console(f"errors: {out_dir / 'errors.log'}")
            tail = _tail_lines(log_file, max_lines=12)
            if tail:
                _console("log tail:")
                for line in tail:
                    _console(line)

        # Mirror per-command progress lines from step logs (quiet mode)
        if not getattr(args, "verbose", False) and not failed:
            if step.name in {"af3score", "flowpacker"}:
                _mirror_log_prefixes(log_file, [f"[{step.name}]"])

        if pending_exc and not getattr(args, "continue_on_error", False):
            raise pending_exc[1].with_traceback(pending_exc[2])

    _console("=" * 78)
    _console("PIPELINE COMPLETE")
    _console("=" * 78)
    if hb:
        hb.complete(extra={"message": "pipeline complete"})
