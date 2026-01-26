from __future__ import annotations

import json
import os
import socket
import threading
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


_HB_CACHE: dict[tuple[str, int], "HeartbeatReporter"] = {}


def _iso(ts: Optional[float]) -> Optional[str]:
    if ts is None:
        return None
    try:
        return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(ts))
    except Exception:
        return None


def _atomic_write_json(path: Path, data: Dict[str, Any]) -> None:
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(path.suffix + ".tmp")
        tmp.write_text(json.dumps(data, ensure_ascii=False, separators=(",", ":")))
        tmp.replace(path)
    except Exception:
        pass


def _get_dist_info() -> Dict[str, int]:
    # Best-effort dist detection
    rank = int(os.environ.get("RANK", "0"))
    world_size = int(os.environ.get("WORLD_SIZE", "1"))
    local_rank = int(os.environ.get("LOCAL_RANK", "0"))
    try:
        import torch

        if torch.distributed.is_available() and torch.distributed.is_initialized():
            rank = int(torch.distributed.get_rank())
            world_size = int(torch.distributed.get_world_size())
            local_rank = int(os.environ.get("LOCAL_RANK", local_rank))
    except Exception:
        pass
    return {"rank": rank, "world_size": world_size, "local_rank": local_rank}


class HeartbeatReporter:
    def __init__(
        self,
        output_dir: str | Path,
        *,
        interval_seconds: Optional[float] = None,
        throughput_window_seconds: float = 60.0,
    ) -> None:
        self.output_dir = Path(output_dir)
        self.interval_seconds = (
            float(os.environ.get("PPIFLOW_HEARTBEAT_INTERVAL", "30"))
            if interval_seconds is None
            else float(interval_seconds)
        )
        if self.interval_seconds <= 0:
            self.interval_seconds = 30.0
        self.throughput_window_seconds = float(throughput_window_seconds)

        d = _get_dist_info()
        self.rank = d["rank"]
        self.world_size = d["world_size"]
        self.local_rank = d["local_rank"]

        self._start_ts: float = 0.0
        self._last_write_ts: float = 0.0
        self._expected_total: Optional[int] = None
        self._produced_total: int = 0
        self._recent: List[Tuple[float, int]] = []
        self._state: str = "running"
        self._primary_counter: str = "outputs"

    @classmethod
    def from_env(cls, output_dir: str | Path) -> Optional["HeartbeatReporter"]:
        if os.environ.get("PPIFLOW_HEARTBEAT", "1").strip() in {"0", "false", "False"}:
            return None
        d = _get_dist_info()
        key = (str(Path(output_dir).expanduser().resolve()), int(d["rank"]))
        hb = _HB_CACHE.get(key)
        if hb is None:
            hb = cls(output_dir)
            _HB_CACHE[key] = hb
        return hb

    def _rank_path(self) -> Path:
        return self.output_dir / f"status_rank{self.rank}.json"

    def _global_path(self) -> Path:
        return self.output_dir / "status.json"

    def start(self, *, expected_total: Optional[int] = None, primary_counter: Optional[str] = None) -> None:
        self._start_ts = time.time()
        self._recent = [(self._start_ts, 0)]
        if expected_total is not None:
            self._expected_total = max(int(expected_total), 0)
        if primary_counter:
            self._primary_counter = str(primary_counter)
        self.update(produced_total=0, expected_total=self._expected_total, force=True)

    def update(
        self,
        *,
        produced_total: int,
        expected_total: Optional[int] = None,
        extra: Optional[Dict[str, Any]] = None,
        state: Optional[str] = None,
        primary_counter: Optional[str] = None,
        counters: Optional[list[dict]] = None,
        last_output_path: Optional[str] = None,
        last_write_time: Optional[float] = None,
        force: bool = False,
    ) -> None:
        now = time.time()
        if not force and (now - self._last_write_ts) < self.interval_seconds:
            return
        self._produced_total = max(int(produced_total), 0)
        if expected_total is not None:
            self._expected_total = max(int(expected_total), 0)
        self._state = state or "running"
        if primary_counter:
            self._primary_counter = str(primary_counter)

        self._recent.append((now, self._produced_total))
        cutoff = now - self.throughput_window_seconds
        while len(self._recent) > 1 and self._recent[0][0] < cutoff:
            self._recent.pop(0)

        if len(self._recent) >= 2:
            dt = max(self._recent[-1][0] - self._recent[0][0], 1e-6)
            dd = max(self._recent[-1][1] - self._recent[0][1], 0)
            rate_per_sec = dd / dt
        else:
            rate_per_sec = 0.0

        expected = max(int(self._expected_total or 0), 1)
        produced = max(int(self._produced_total), 0)
        percent = min(max(produced / expected, 0.0), 1.0)

        remaining = max(expected - produced, 0)
        eta_seconds: Optional[float] = None
        if state == "completed":
            percent = 1.0
            eta_seconds = 0.0
        elif rate_per_sec > 0:
            eta_seconds = remaining / rate_per_sec
        else:
            elapsed = max(now - (self._start_ts or now), 1e-6)
            overall_rate = produced / elapsed
            if overall_rate > 0:
                eta_seconds = remaining / overall_rate

        stage = os.environ.get("PPIFLOW_STAGE", None)
        step = os.environ.get("PPIFLOW_STEP", None)
        protocol = os.environ.get("PPIFLOW_PROTOCOL", None)
        run_id = os.environ.get("PPIFLOW_RUN_ID", None)
        task_name = os.environ.get("PPIFLOW_TASK_NAME", None)
        seed = os.environ.get("PPIFLOW_SEED", None)

        payload: Dict[str, Any] = {
            "job": {
                "output_dir": str(self.output_dir.resolve()),
                "pid": os.getpid(),
                "host": socket.gethostname(),
            },
            "pipeline": {
                "protocol": protocol,
                "run_id": int(run_id) if str(run_id).isdigit() else run_id,
                "stage": stage,
                "step": step,
                "task_name": task_name,
                "seed": int(seed) if str(seed).lstrip("-").isdigit() else seed,
            },
            "status": {
                "state": self._state,
                "started_at": _iso(self._start_ts or None),
                "updated_at": _iso(now),
                "percent": percent,
                "eta_seconds": eta_seconds,
                "eta_timestamp": _iso(now + eta_seconds) if isinstance(eta_seconds, (int, float)) else None,
            },
            "compute": {
                "global_rank": self.rank,
                "local_rank": self.local_rank,
                "world_size": self.world_size,
            },
            "progress": {
                "expected_total": expected,
                "produced_total": produced,
                "produced_since_start": max(produced - int(self._recent[0][1]), 0) if self._recent else produced,
                "throughput_per_min": rate_per_sec * 60.0,
                "throughput_window_sec": self.throughput_window_seconds,
                "primary_counter": self._primary_counter,
                "counters": counters or [],
            },
            "writer": {
                "last_output_path": last_output_path,
                "last_write_time": _iso(last_write_time) if last_write_time else None,
            },
        }
        if extra is not None:
            payload["extra"] = extra

        _atomic_write_json(self._rank_path(), payload)

        if self.rank == 0:
            agg = self._aggregate(now=now)
            _atomic_write_json(self._global_path(), agg)

        self._last_write_ts = now

    def touch(self, *, extra: Optional[Dict[str, Any]] = None, force: bool = False) -> None:
        self.update(
            produced_total=self._produced_total,
            expected_total=self._expected_total,
            extra=extra,
            state=self._state,
            force=force,
        )

    def complete(self, *, extra: Optional[Dict[str, Any]] = None) -> None:
        self.update(
            produced_total=self._produced_total,
            expected_total=self._expected_total,
            extra=extra,
            state="completed",
            force=True,
        )

    def _aggregate(self, *, now: float) -> Dict[str, Any]:
        per_rank: List[Dict[str, Any]] = []
        produced_sum = 0
        expected_sum = 0

        for r in range(max(int(self.world_size), 1)):
            p = self.output_dir / f"status_rank{r}.json"
            try:
                d = json.loads(p.read_text())
                prog = d.get("progress", {}) if isinstance(d, dict) else {}
                produced = int(prog.get("produced_total", 0) or 0)
                expected = int(prog.get("expected_total", 0) or 0)
                produced_sum += max(produced, 0)
                expected_sum += max(expected, 0)
                per_rank.append(
                    {
                        "rank": r,
                        "produced_total": produced,
                        "expected_total": expected,
                        "updated_at": (d.get("status", {}) or {}).get("updated_at"),
                    }
                )
            except Exception:
                per_rank.append(
                    {"rank": r, "produced_total": 0, "expected_total": 0, "updated_at": None}
                )

        expected = max(expected_sum, 1)
        percent = min(max(produced_sum / expected, 0.0), 1.0)

        stage = os.environ.get("PPIFLOW_STAGE", None)
        step = os.environ.get("PPIFLOW_STEP", None)
        protocol = os.environ.get("PPIFLOW_PROTOCOL", None)
        run_id = os.environ.get("PPIFLOW_RUN_ID", None)
        task_name = os.environ.get("PPIFLOW_TASK_NAME", None)
        seed = os.environ.get("PPIFLOW_SEED", None)

        return {
            "job": {
                "output_dir": str(self.output_dir.resolve()),
                "pid": os.getpid(),
                "host": socket.gethostname(),
            },
            "pipeline": {
                "protocol": protocol,
                "run_id": int(run_id) if str(run_id).isdigit() else run_id,
                "stage": stage,
                "step": step,
                "task_name": task_name,
                "seed": int(seed) if str(seed).lstrip("-").isdigit() else seed,
            },
            "status": {
                "state": self._state,
                "started_at": _iso(self._start_ts or None),
                "updated_at": _iso(now),
                "percent": percent,
            },
            "compute": {
                "global_rank": 0,
                "world_size": self.world_size,
            },
            "progress": {
                "expected_total": expected_sum,
                "produced_total": produced_sum,
                "primary_counter": self._primary_counter,
            },
            "distributed": {"per_rank": per_rank},
        }


def start_keepalive(
    hb: Optional[HeartbeatReporter],
    *,
    interval_s: Optional[float] = None,
    extra: Optional[Dict[str, Any]] = None,
) -> Optional[tuple[threading.Event, threading.Thread]]:
    if hb is None:
        return None
    if interval_s is None:
        interval_s = float(os.environ.get("PPIFLOW_STEP_HEARTBEAT_INTERVAL", "30") or 30)
    if interval_s <= 0:
        return None

    stop = threading.Event()

    def _loop():
        while not stop.wait(interval_s):
            try:
                hb.touch(extra=extra)
            except Exception:
                pass

    thread = threading.Thread(target=_loop, daemon=True)
    thread.start()
    return stop, thread
