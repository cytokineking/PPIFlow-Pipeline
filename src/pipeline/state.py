from __future__ import annotations

import hashlib
import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional

from . import PIPELINE_VERSION
from .io import read_json, write_json


class StateError(RuntimeError):
    pass


def sha256_json(data: Any) -> str:
    payload = json.dumps(data, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _tool_stamp(path: Optional[str]) -> Optional[dict]:
    if not path:
        return None
    p = Path(path)
    if not p.exists():
        return {"path": str(p), "exists": False}
    stamp = {"path": str(p), "exists": True}
    if p.is_file():
        stamp["sha256"] = sha256_file(p)
        stamp["mtime"] = int(p.stat().st_mtime)
        stamp["size"] = int(p.stat().st_size)
    else:
        stamp["mtime"] = int(p.stat().st_mtime)
    return stamp


def collect_tool_versions(tools: dict) -> dict:
    versions = {}
    for key, path in (tools or {}).items():
        versions[key] = _tool_stamp(path)
    return versions


def load_state(path: str | Path) -> Optional[dict]:
    try:
        return read_json(path)
    except Exception:
        return None


def save_state(path: str | Path, data: dict) -> None:
    write_json(path, data, indent=2)


def validate_state(state: dict, input_sha256: str, tool_versions: dict) -> None:
    job = state.get("job") or {}
    if job.get("input_sha256") and job.get("input_sha256") != input_sha256:
        raise StateError(
            "pipeline_state.json input hash does not match pipeline_input.json. "
            "Use a new output directory."
        )
    prev_tools = job.get("tool_versions") or {}
    for key, prev in prev_tools.items():
        cur = tool_versions.get(key)
        if prev is None or cur is None:
            continue
        if prev.get("path") != cur.get("path"):
            raise StateError(
                f"pipeline_state.json tool path mismatch for {key}. Use a new output directory."
            )
        if prev.get("sha256") and cur.get("sha256") and prev.get("sha256") != cur.get("sha256"):
            raise StateError(
                f"pipeline_state.json tool checksum mismatch for {key}. Use a new output directory."
            )


def rebuild_state(out_dir: str | Path, *, input_sha256: str, tool_versions: dict, target_n: int) -> dict:
    runs = [{"run_id": 0, "run_seed": int(time.time_ns() % (2**31 - 1)), "target_samples": int(target_n)}]
    return {
        "version": PIPELINE_VERSION,
        "layout": "v1",
        "job": {
            "created_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "input_sha256": input_sha256,
            "pipeline_version": PIPELINE_VERSION,
            "tool_versions": tool_versions,
            "recovered": True,
        },
        "runs": runs,
    }


def init_or_update_state(
    *,
    out_dir: str | Path,
    input_sha256: str,
    tool_versions: dict,
    target_n: int,
    seeds: Optional[list[int]] = None,
) -> dict:
    state_path = Path(out_dir) / "pipeline_state.json"
    state = load_state(state_path)
    if state is None:
        state = rebuild_state(out_dir, input_sha256=input_sha256, tool_versions=tool_versions, target_n=target_n)
    else:
        if state.get("version") != PIPELINE_VERSION:
            raise StateError("pipeline_state.json version mismatch. Use a new output directory.")
        validate_state(state, input_sha256, tool_versions)

    runs = list(state.get("runs") or [])
    if not runs:
        runs = [{"run_id": 0, "run_seed": int(time.time_ns() % (2**31 - 1)), "target_samples": int(target_n)}]
    # Ensure run 0 exists
    if runs and runs[0].get("run_id") != 0:
        runs.insert(0, {"run_id": 0, "run_seed": int(time.time_ns() % (2**31 - 1)), "target_samples": int(target_n)})
    if seeds:
        if len(seeds) < 1:
            raise StateError("seeds list is empty")
        runs[0]["run_seed"] = int(seeds[0])
    # Update tool versions
    state["job"] = state.get("job") or {}
    state["job"].update({
        "input_sha256": input_sha256,
        "pipeline_version": PIPELINE_VERSION,
        "tool_versions": tool_versions,
    })
    state["runs"] = runs
    save_state(state_path, state)
    return state
