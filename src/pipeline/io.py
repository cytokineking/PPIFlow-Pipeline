from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import yaml


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def _atomic_write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(content)
    tmp.replace(path)


def write_json(path: str | Path, data: Any, *, indent: int | None = 2) -> None:
    _atomic_write(Path(path), json.dumps(data, indent=indent, sort_keys=False))


def read_json(path: str | Path) -> Any:
    with open(path, "r") as f:
        return json.load(f)


def write_yaml(path: str | Path, data: Any) -> None:
    _atomic_write(Path(path), yaml.dump(data, sort_keys=False, default_flow_style=False))


def read_yaml(path: str | Path) -> Any:
    with open(path, "r") as f:
        return yaml.safe_load(f)


def resolve_path(path: str | Path, *, base_dir: str | Path | None = None) -> str:
    p = Path(path)
    if not p.is_absolute():
        if base_dir is None:
            p = (Path.cwd() / p).resolve()
        else:
            p = (Path(base_dir) / p).resolve()
    else:
        p = p.resolve()
    return str(p)


def resolve_optional_path(path: str | None, *, base_dir: str | Path | None = None) -> str | None:
    if not path:
        return None
    return resolve_path(path, base_dir=base_dir)


def list_files(path: str | Path, pattern: str) -> list[Path]:
    p = Path(path)
    if not p.exists():
        return []
    return list(p.glob(pattern))


def safe_relpath(path: str | Path, base: str | Path) -> str:
    try:
        return str(Path(path).resolve().relative_to(Path(base).resolve()))
    except Exception:
        return str(path)


def repo_root(start: str | Path | None = None) -> Path:
    if start is None:
        start_path = Path(__file__).resolve()
    else:
        start_path = Path(start).resolve()
    for parent in [start_path] + list(start_path.parents):
        if (parent / "ppiflow.py").exists() or (parent / "README.md").exists():
            return parent
    return Path.cwd().resolve()
