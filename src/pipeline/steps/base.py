from __future__ import annotations

import csv
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Optional

from ..heartbeat import HeartbeatReporter, start_keepalive
from ..io import ensure_dir


@dataclass
class StepContext:
    out_dir: Path
    input_data: dict
    state: dict
    run_id: int
    rank: int
    world_size: int
    local_rank: int
    reuse: bool
    heartbeat: Optional[HeartbeatReporter]


class StepError(RuntimeError):
    pass


class Step:
    name: str = ""
    stage: str = ""
    supports_indices: bool = True

    def __init__(self, cfg: dict) -> None:
        self.cfg = cfg

    def expected_total(self, ctx: StepContext) -> int:
        raise NotImplementedError

    def scan_done(self, ctx: StepContext) -> set[int]:
        return set()

    def run_indices(self, ctx: StepContext, indices: list[int]) -> None:
        raise NotImplementedError

    def run_full(self, ctx: StepContext) -> None:
        # Default to index-based run
        indices = list(range(self.expected_total(ctx)))
        self.run_indices(ctx, indices)

    def write_manifest(self, ctx: StepContext) -> None:
        return

    def run(self, ctx: StepContext) -> None:
        expected_total = self.expected_total(ctx)
        done = self.scan_done(ctx)
        missing = [i for i in range(expected_total) if i not in done]
        done_count = len(done)
        missing_count = len(missing)
        output_rows = None
        manifest = self.cfg.get("manifest")
        if manifest:
            p = Path(manifest)
            if not p.is_absolute():
                p = ctx.out_dir / p
            if p.exists():
                try:
                    with p.open("r", newline="") as handle:
                        rows = list(csv.reader(handle))
                    output_rows = max(len(rows) - 1, 0)
                except Exception:
                    output_rows = None
        extra = ""
        if output_rows is not None:
            extra = f" outputs={output_rows}"

        if not missing:
            print(
                f"[{self.name}] reuse_check expected={expected_total} done={done_count} missing=0 run=0{extra}",
                file=sys.__stdout__,
                flush=True,
            )
            try:
                if self.cfg.get("manifest"):
                    self.write_manifest(ctx)
            except Exception:
                pass
            return

        if self.supports_indices:
            my_missing = [i for i in missing if (i % max(ctx.world_size, 1)) == ctx.rank]
            if not my_missing:
                print(
                    f"[{self.name}] reuse_check expected={expected_total} done={done_count} missing={missing_count} run=0{extra}",
                    file=sys.__stdout__,
                    flush=True,
                )
                return
            print(
                f"[{self.name}] reuse_check expected={expected_total} done={done_count} missing={missing_count} run={len(my_missing)}{extra}",
                file=sys.__stdout__,
                flush=True,
            )
            owned_total = ((expected_total - 1 - ctx.rank) // max(ctx.world_size, 1) + 1) if expected_total - 1 >= ctx.rank else 0
            hb = ctx.heartbeat
            if hb:
                hb.start(expected_total=owned_total, primary_counter=self.name)
            keepalive = start_keepalive(hb, extra={"step": self.name})
            try:
                self.run_indices(ctx, my_missing)
            finally:
                if keepalive:
                    stop, thread = keepalive
                    stop.set()
                    thread.join(timeout=1.0)
                if hb:
                    hb.complete(extra={"step": self.name})
            try:
                if self.cfg.get("manifest"):
                    self.write_manifest(ctx)
            except Exception:
                pass
        else:
            # Step does not support partial indices; run full step if anything missing.
            print(
                f"[{self.name}] reuse_check expected={expected_total} done={done_count} missing={missing_count} run=1{extra}",
                file=sys.__stdout__,
                flush=True,
            )
            hb = ctx.heartbeat
            if hb:
                hb.start(expected_total=expected_total, primary_counter=self.name)
            keepalive = start_keepalive(hb, extra={"step": self.name})
            try:
                self.run_full(ctx)
            finally:
                if keepalive:
                    stop, thread = keepalive
                    stop.set()
                    thread.join(timeout=1.0)
                if hb:
                    hb.complete(extra={"step": self.name})
            try:
                if self.cfg.get("manifest"):
                    self.write_manifest(ctx)
            except Exception:
                pass

    def output_dir(self, ctx: StepContext) -> Path:
        out_dir = self.cfg.get("output_dir")
        if not out_dir:
            raise StepError(f"Step {self.name} missing output_dir in config")
        p = Path(out_dir)
        if not p.is_absolute():
            p = ctx.out_dir / p
        ensure_dir(p)
        return p

    def manifest_path(self, ctx: StepContext) -> Path:
        manifest = self.cfg.get("manifest")
        if not manifest:
            raise StepError(f"Step {self.name} missing manifest path in config")
        p = Path(manifest)
        if not p.is_absolute():
            p = ctx.out_dir / p
        ensure_dir(p.parent)
        return p
