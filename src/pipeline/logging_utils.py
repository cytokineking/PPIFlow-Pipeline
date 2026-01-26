from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path
from typing import Iterable


def run_command(
    cmd: Iterable[str],
    *,
    env: dict | None = None,
    cwd: str | Path | None = None,
    log_file: str | Path | None = None,
    verbose: bool = False,
) -> None:
    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        if verbose:
            with log_path.open("a") as handle:
                proc = subprocess.Popen(
                    list(cmd),
                    cwd=str(cwd) if cwd else None,
                    env=env,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1,
                )
                assert proc.stdout is not None
                for line in proc.stdout:
                    sys.stdout.write(line)
                    handle.write(line)
                ret = proc.wait()
                if ret != 0:
                    raise subprocess.CalledProcessError(ret, list(cmd))
        else:
            with log_path.open("a") as handle:
                subprocess.check_call(
                    list(cmd),
                    cwd=str(cwd) if cwd else None,
                    env=env,
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                )
    else:
        subprocess.check_call(list(cmd), cwd=str(cwd) if cwd else None, env=env)


def log_command_progress(
    step: str,
    idx: int,
    total: int,
    *,
    item: str | None = None,
    phase: str | None = None,
    status: str = "OK",
    elapsed: float | None = None,
    log_file: str | Path | None = None,
    extra: str | None = None,
) -> None:
    parts = [f"[{step}] {idx}/{total}", status]
    if elapsed is not None:
        parts.append(f"elapsed={elapsed:.2f}s")
    if item:
        parts.append(f"item={item}")
    if phase:
        parts.append(f"phase={phase}")
    if extra:
        parts.append(extra)
    if status != "OK" and log_file:
        parts.append(f"log={log_file}")
    print(" ".join(parts), file=sys.__stdout__, flush=True)
