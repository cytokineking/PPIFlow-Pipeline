#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

from pipeline.hotspots import expand_hotspots, parse_chain_list


def main() -> None:
    pdb_path = Path(__file__).resolve().parents[1] / "assets" / "frameworks" / "7xl0_nanobody_framework.pdb"
    chains = parse_chain_list("A,C")
    assert chains == ["A", "C"], chains

    expanded = expand_hotspots("A", pdb_path=pdb_path)
    assert expanded, "Chain-only expansion returned empty list"
    assert all(t.startswith("A") for t in expanded), expanded[:5]

    expanded_range = expand_hotspots("A1-3", pdb_path=pdb_path)
    assert expanded_range[:3] == ["A1", "A2", "A3"], expanded_range

    expanded_mixed = expand_hotspots(["A1-2", "A5"], pdb_path=pdb_path)
    assert expanded_mixed == ["A1", "A2", "A5"], expanded_mixed

    print("hotspot parser tests passed")


if __name__ == "__main__":
    main()
