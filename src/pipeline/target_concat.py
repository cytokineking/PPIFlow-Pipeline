from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

from Bio import PDB

from .io import ensure_dir, write_json


_HOTSPOT_RE = re.compile(r"^([A-Za-z])(\d+)([A-Za-z]?)$")
_MAX_PDB_RESSEQ = 9999
_MAX_HOTSPOT_CLI_LEN = 8000


def concatenate_target_chains(
    pdb_path: str | Path,
    target_chains: Iterable[str],
    output_dir: str | Path,
    *,
    gap_residues: int = 50,
) -> dict:
    pdb_path = Path(pdb_path)
    out_root = Path(output_dir) / "inputs"
    ensure_dir(out_root)

    if gap_residues < 0:
        raise ValueError("gap_residues must be non-negative")

    concat_path = out_root / "concatenated_target.pdb"
    chain_map_path = out_root / "target_chain_map.json"
    offsets_path = out_root / "target_chain_offsets.json"

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("target", str(pdb_path))
    model = structure[0]

    chain_order = [str(c) for c in target_chains]
    if not chain_order:
        raise ValueError("target_chains is empty")

    new_structure = PDB.Structure.Structure("target_concat")
    new_model = PDB.Model.Model(0)
    new_chain = PDB.Chain.Chain("B")

    mapping_entries: list[dict[str, object]] = []
    offsets: list[dict[str, object]] = []

    next_resseq = 1
    for idx, chain_id in enumerate(chain_order):
        if not model.has_id(chain_id):
            raise ValueError(f"Target chain not found in PDB: {chain_id}")
        chain = model[chain_id]
        start_resseq = next_resseq
        count = 0
        for res in chain:
            if not PDB.is_aa(res, standard=False):
                continue
            res_id = res.get_id()
            orig_resseq = int(res_id[1])
            orig_icode = str(res_id[2]).strip()
            new_resseq = next_resseq
            if new_resseq > _MAX_PDB_RESSEQ:
                raise ValueError(
                    f"Concatenated target exceeds PDB resseq limit ({_MAX_PDB_RESSEQ})."
                )
            next_resseq += 1
            count += 1

            new_res = PDB.Residue.Residue((" ", new_resseq, " "), res.resname, res.segid)
            for atom in res:
                new_res.add(atom.copy())
            new_chain.add(new_res)

            mapping_entries.append({
                "chain": chain_id,
                "resseq": orig_resseq,
                "icode": orig_icode,
                "resseq_B": new_resseq,
            })

        end_resseq = next_resseq - 1
        offsets.append({
            "chain": chain_id,
            "offset": start_resseq - 1,
            "length": count,
            "start_resseq_B": start_resseq,
            "end_resseq_B": end_resseq,
        })
        if idx < len(chain_order) - 1 and count > 0:
            next_resseq += int(gap_residues)

    if not mapping_entries:
        raise ValueError("No residues copied while concatenating target chains.")

    new_model.add(new_chain)
    new_structure.add(new_model)
    io = PDB.PDBIO()
    io.set_structure(new_structure)
    io.save(str(concat_path))

    map_payload = {
        "entries": mapping_entries,
        "gap_residues": int(gap_residues),
        "target_chain_order": chain_order,
        "source_pdb": str(pdb_path),
    }
    offsets_payload = {
        "chains": offsets,
        "gap_residues": int(gap_residues),
        "target_chain_order": chain_order,
        "source_pdb": str(pdb_path),
    }
    write_json(chain_map_path, map_payload, indent=2)
    write_json(offsets_path, offsets_payload, indent=2)

    return {
        "concatenated_pdb": str(concat_path),
        "chain_map_path": str(chain_map_path),
        "offsets_path": str(offsets_path),
        "chain_map_entries": mapping_entries,
    }


def map_hotspots_to_concatenated(
    hotspots: Iterable[str],
    chain_map_entries: Iterable[dict],
) -> list[str]:
    mapping: dict[tuple[str, int, str], int] = {}
    for entry in chain_map_entries:
        chain = str(entry.get("chain") or "")
        resseq = int(entry.get("resseq") or 0)
        icode = str(entry.get("icode") or "").strip()
        resseq_b = int(entry.get("resseq_B") or 0)
        mapping[(chain, resseq, icode)] = resseq_b

    mapped: list[int] = []
    for token in hotspots:
        token = str(token).strip()
        if not token:
            continue
        match = _HOTSPOT_RE.match(token)
        if not match:
            raise ValueError(f"Unrecognized hotspot token: {token}")
        chain_id = match.group(1)
        resseq = int(match.group(2))
        icode = match.group(3) or ""
        key = (chain_id, resseq, icode)
        if key not in mapping:
            raise ValueError(f"Hotspot not found in chain map: {token}")
        mapped.append(mapping[key])

    mapped = sorted(set(mapped))
    return [f"B{r}" for r in mapped]


def compress_hotspots(hotspots: Iterable[str]) -> str:
    resseqs: list[int] = []
    for token in hotspots:
        token = str(token).strip()
        if not token:
            continue
        if not token.startswith("B"):
            raise ValueError(f"Expected chain B hotspot token, got: {token}")
        resseqs.append(int(token[1:]))
    if not resseqs:
        return ""
    resseqs = sorted(set(resseqs))
    ranges: list[tuple[int, int]] = []
    start = prev = resseqs[0]
    for r in resseqs[1:]:
        if r == prev + 1:
            prev = r
            continue
        ranges.append((start, prev))
        start = prev = r
    ranges.append((start, prev))
    parts: list[str] = []
    for a, b in ranges:
        if a == b:
            parts.append(f"B{a}")
        else:
            parts.append(f"B{a}-B{b}")
    return ",".join(parts)


def maybe_write_hotspots_file(
    compressed_hotspots: str,
    output_dir: str | Path,
) -> tuple[str | None, str | None]:
    if not compressed_hotspots:
        return None, None
    if len(compressed_hotspots) <= _MAX_HOTSPOT_CLI_LEN:
        return compressed_hotspots, None
    out_root = Path(output_dir) / "inputs"
    ensure_dir(out_root)
    path = out_root / "hotspots_mapped.txt"
    path.write_text(compressed_hotspots + "\n")
    return None, str(path)
