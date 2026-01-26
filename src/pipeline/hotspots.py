from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable, List

from Bio import PDB


_CHAIN_RE = re.compile(r"^[A-Za-z]$")
_SINGLE_RE = re.compile(r"^([A-Za-z])(\d+)([A-Za-z]?)$")
_RANGE_RE = re.compile(r"^([A-Za-z])(\d+)-([A-Za-z])?(\d+)$")


def parse_chain_list(value: object | None) -> List[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        out: List[str] = []
        for item in value:
            out.extend(parse_chain_list(item))
        return out
    text = str(value).strip()
    if not text:
        return []
    if "," in text or "_" in text:
        parts = re.split(r"[,_]", text)
        return [p.strip() for p in parts if p.strip()]
    # Treat concatenated chains as a list of single-letter chain IDs.
    return [c for c in text if c.strip()]


def _tokenize_hotspots(value: object | None) -> List[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        tokens: List[str] = []
        for item in value:
            tokens.extend(_tokenize_hotspots(item))
        return tokens
    text = str(value).strip()
    if not text:
        return []
    return [t.strip() for t in text.split(",") if t.strip()]


def _read_hotspots_file(path: str | Path) -> List[str]:
    p = Path(path)
    text = p.read_text().strip()
    if not text:
        return []
    # Allow JSON list syntax as a convenience.
    if text.lstrip().startswith("[") and text.rstrip().endswith("]"):
        try:
            import json

            data = json.loads(text)
            if isinstance(data, list):
                return [str(x).strip() for x in data if str(x).strip()]
        except Exception:
            pass
    tokens = re.split(r"[,\s]+", text)
    return [t.strip() for t in tokens if t.strip()]


def resolve_hotspots_input(
    value: object | None,
    *,
    hotspots_file: str | Path | None = None,
) -> object | None:
    if hotspots_file:
        return _read_hotspots_file(hotspots_file)
    if isinstance(value, str):
        text = value.strip()
        if text.startswith("@") and len(text) > 1:
            return _read_hotspots_file(text[1:])
    return value


def _chain_residue_ids_from_pdb(pdb_path: str | Path, chain_id: str) -> List[str]:
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("hotspots", str(pdb_path))
    model = structure[0]
    if not model.has_id(chain_id):
        return []
    chain = model[chain_id]
    ids: List[str] = []
    for res in chain:
        if not PDB.is_aa(res, standard=True):
            continue
        res_num = str(res.get_id()[1])
        icode = res.get_id()[2].strip()
        ids.append(f"{chain_id}{res_num}{icode}")
    return ids


def expand_hotspots(
    value: object | None,
    *,
    pdb_path: str | Path | None = None,
) -> List[str]:
    tokens = _tokenize_hotspots(value)
    if not tokens:
        return []

    expanded: List[str] = []
    chain_only: List[str] = []

    for token in tokens:
        if _CHAIN_RE.match(token):
            chain_only.append(token)
            continue

        m_range = _RANGE_RE.match(token)
        if m_range:
            chain_a = m_range.group(1)
            start = int(m_range.group(2))
            chain_b = m_range.group(3)
            end = int(m_range.group(4))
            if chain_b and chain_b != chain_a:
                raise ValueError(f"Hotspot range spans multiple chains: {token}")
            if end < start:
                raise ValueError(f"Hotspot range must be ascending: {token}")
            expanded.extend([f"{chain_a}{i}" for i in range(start, end + 1)])
            continue

        m_single = _SINGLE_RE.match(token)
        if m_single:
            expanded.append(token)
            continue

        raise ValueError(f"Unrecognized hotspot token: {token}")

    if chain_only:
        if pdb_path is None:
            raise ValueError("Chain-only hotspots require pdb_path to expand")
        for chain_id in chain_only:
            expanded.extend(_chain_residue_ids_from_pdb(pdb_path, chain_id))

    return expanded
