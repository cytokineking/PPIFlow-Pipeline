from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

from Bio import PDB

from .config import InputSpec, apply_preset, build_input_from_cli, load_input, normalize_input, validate_input, write_cli_input_yaml
from .io import ensure_dir, write_json, write_yaml, repo_root
from .state import collect_tool_versions, init_or_update_state, sha256_json
from .steps import STEP_ORDER


def _step_config(name: str, run_id: int, input_data: dict) -> dict:
    run_dir = "output"
    manifests_dir = f"{run_dir}/manifests"
    cfg: dict[str, Any] = {
        "name": name,
        "stage": None,
        "output_dir": None,
        "manifest": None,
        "command": None,
        "input_dir": None,
    }
    if name == "gen":
        cfg.update({
            "stage": "gen",
            "output_dir": f"{run_dir}/backbones",
            "manifest": f"{manifests_dir}/backbones.csv",
        })
    elif name == "seq1":
        cfg.update({
            "stage": "seq",
            "output_dir": f"{run_dir}/seqs_round1",
            "manifest": f"{manifests_dir}/seqs_round1.csv",
            "input_dir": f"{run_dir}/backbones",
        })
    elif name == "flowpacker1":
        cfg.update({
            "stage": "score",
            "output_dir": f"{run_dir}/flowpacker_round1",
            "manifest": f"{manifests_dir}/flowpacker_round1.csv",
            "input_dir": f"{run_dir}/seqs_round1",
        })
    elif name == "af3score1":
        cfg.update({
            "stage": "score",
            "output_dir": f"{run_dir}/af3score_round1",
            "manifest": f"{manifests_dir}/af3score_round1.csv",
            "input_dir": f"{run_dir}/flowpacker_round1",
        })
    elif name == "rosetta_interface":
        af3_r1 = ((input_data.get("filters") or {}).get("af3score") or {}).get("round1") or {}
        use_filtered = af3_r1.get("iptm_min") is not None or af3_r1.get("top_k") is not None
        cfg.update({
            "stage": "rosetta",
            "output_dir": f"{run_dir}/rosetta_interface",
            "manifest": f"{manifests_dir}/rosetta_interface.csv",
            "input_dir": f"{run_dir}/af3score_round1/filtered_pdbs" if use_filtered else f"{run_dir}/af3score_round1",
        })
    elif name == "interface_enrich":
        cfg.update({
            "stage": "rosetta",
            "output_dir": f"{run_dir}/interface_enrich",
            "manifest": f"{run_dir}/interface_enrich/fixed_positions.csv",
            "residue_energy_csv": f"{run_dir}/rosetta_interface/residue_energy.csv",
        })
    elif name == "partial":
        cfg.update({
            "stage": "partial",
            "output_dir": f"{run_dir}/partial_flow",
            "manifest": f"{manifests_dir}/partial_flow.csv",
            "fixed_positions_csv": f"{run_dir}/interface_enrich/fixed_positions.csv",
        })
    elif name == "seq2":
        cfg.update({
            "stage": "seq",
            "output_dir": f"{run_dir}/seqs_round2",
            "manifest": f"{manifests_dir}/seqs_round2.csv",
            "input_dir": f"{run_dir}/partial_flow",
            "fixed_positions_csv": f"{run_dir}/interface_enrich/fixed_positions.csv",
        })
    elif name == "flowpacker2":
        cfg.update({
            "stage": "score",
            "output_dir": f"{run_dir}/flowpacker_round2",
            "manifest": f"{manifests_dir}/flowpacker_round2.csv",
            "input_dir": f"{run_dir}/seqs_round2",
        })
    elif name == "af3score2":
        cfg.update({
            "stage": "score",
            "output_dir": f"{run_dir}/af3score_round2",
            "manifest": f"{manifests_dir}/af3score_round2.csv",
            "input_dir": f"{run_dir}/flowpacker_round2",
        })
    elif name == "relax":
        af3_r2 = ((input_data.get("filters") or {}).get("af3score") or {}).get("round2") or {}
        use_filtered = af3_r2.get("iptm_min") is not None or af3_r2.get("ptm_min") is not None or af3_r2.get("top_k") is not None
        cfg.update({
            "stage": "rosetta",
            "output_dir": f"{run_dir}/relax",
            "manifest": f"{manifests_dir}/relax.csv",
            "input_dir": f"{run_dir}/af3score_round2/filtered_pdbs" if use_filtered else f"{run_dir}/af3score_round2",
        })
    elif name == "dockq":
        tools = input_data.get("tools") or {}
        use_refold = bool(
            input_data.get("optional_af3_refold")
            and tools.get("af3score_repo")
            and tools.get("af3_weights")
        )
        cfg.update({
            "stage": "score",
            "output_dir": f"{run_dir}/dockq",
            "manifest": f"{manifests_dir}/dockq.csv",
            "input_dir": f"{run_dir}/af3_refold/pdbs" if use_refold else f"{run_dir}/relax",
        })
    elif name == "rank":
        cfg.update({
            "stage": "rank",
            "output_dir": "results",
            "manifest": "results/summary.csv",
        })
        cfg["ranking"] = input_data.get("ranking") or {}
    elif name == "af3_refold":
        cfg.update({
            "stage": "score",
            "output_dir": f"{run_dir}/af3_refold",
            "manifest": f"{manifests_dir}/af3_refold.csv",
            "input_dir": f"{run_dir}/relax",
        })
    return cfg


def _available_chain_ids() -> list[str]:
    return list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")


def _rename_framework_chains(pdb_path: Path, mapping: dict[str, str], out_path: Path) -> None:
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("framework", str(pdb_path))
    for model in structure:
        for chain in model:
            if chain.id in mapping:
                chain.id = mapping[chain.id]
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(str(out_path))


def _auto_rename_framework(input_data: dict, out_dir: Path) -> dict:
    protocol = input_data.get("protocol")
    if protocol == "binder":
        out = dict(input_data)
        out["binder_chain"] = "A"
        return out
    if protocol not in {"antibody", "vhh"}:
        return input_data
    framework = input_data.get("framework") or {}
    heavy = framework.get("heavy_chain")
    light = framework.get("light_chain")
    if not heavy:
        return input_data

    framework_pdb = framework.get("pdb")
    if not framework_pdb:
        return input_data

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("framework", str(framework_pdb))
    model = structure[0]
    chain_ids = [c.id for c in model.get_chains()]

    desired: dict[str, str] = {}
    for cid in chain_ids:
        if cid == heavy:
            desired[cid] = "A"
        elif light and cid == light:
            desired[cid] = "C"
        else:
            desired[cid] = cid

    reserved = {"A", "C", "B"}
    mapping: dict[str, str] = {}
    used: set[str] = set()

    # Reserve canonical IDs for heavy/light first.
    for orig, target in desired.items():
        if target in {"A", "C"}:
            mapping[orig] = target
            used.add(target)

    for orig, target in desired.items():
        if orig in mapping:
            continue
        if target not in used and target not in reserved:
            mapping[orig] = target
            used.add(target)
            continue
        for cid in _available_chain_ids():
            if cid in used or cid in reserved:
                continue
            mapping[orig] = cid
            used.add(cid)
            break

    if not mapping or all(orig == new for orig, new in mapping.items()):
        return input_data

    out_inputs = out_dir / "inputs"
    ensure_dir(out_inputs)
    renamed_path = out_inputs / f"framework_renamed_{Path(framework_pdb).stem}.pdb"
    _rename_framework_chains(Path(framework_pdb), mapping, renamed_path)

    framework = dict(framework)
    framework["pdb"] = str(renamed_path)
    framework["heavy_chain"] = mapping.get(heavy, "A")
    if light:
        framework["light_chain"] = mapping.get(light, "C")
    out = dict(input_data)
    out["framework"] = framework
    return out


def _resolve_binder_chain(input_data: dict) -> str:
    protocol = input_data.get("protocol")
    if protocol in {"antibody", "vhh"}:
        framework = input_data.get("framework") or {}
        return str(framework.get("heavy_chain") or "A")
    return str(input_data.get("binder_chain") or "A")


def _apply_default_command(step_name: str, run_id: int, out_dir: Path, input_data: dict, cfg: dict) -> None:
    tools = input_data.get("tools") or {}
    root = repo_root()

    if step_name in {"flowpacker1", "flowpacker2"}:
        flowpacker_repo = tools.get("flowpacker_repo")
        if not flowpacker_repo:
            return
        script = root / "scripts" / "run_flowpacker.py"
        if not script.exists():
            return
        run_dir = out_dir / "output"
        binder_chain = _resolve_binder_chain(input_data)
        num_jobs = int(os.environ.get("PPIFLOW_FLOWPACKER_JOBS", "1"))
        if step_name == "flowpacker1":
            input_pdb_dir = run_dir / "backbones"
            seq_fasta_dir = run_dir / "seqs_round1" / "seqs"
            output_dir = run_dir / "flowpacker_round1"
        else:
            input_pdb_dir = run_dir / "seqs_round2" / "pdbs"
            seq_fasta_dir = run_dir / "seqs_round2" / "seqs"
            output_dir = run_dir / "flowpacker_round2"
        cfg["command"] = [
            "python",
            str(script),
            "--input_pdb_dir",
            str(input_pdb_dir),
            "--seq_fasta_dir",
            str(seq_fasta_dir),
            "--output_dir",
            str(output_dir),
            "--flowpacker_repo",
            str(flowpacker_repo),
            "--binder_chain",
            str(binder_chain),
            "--num_jobs",
            str(num_jobs),
        ]
        return

    if step_name in {"af3score1", "af3score2"}:
        af3_repo = tools.get("af3score_repo")
        af3_weights = tools.get("af3_weights")
        if not af3_repo or not af3_weights:
            return
        script = root / "scripts" / "run_af3score.py"
        if not script.exists():
            return
        run_dir = out_dir / "output"
        num_jobs = int(os.environ.get("PPIFLOW_AF3_JOBS", "1"))
        if step_name == "af3score1":
            input_pdb_dir = run_dir / "flowpacker_round1" / "flowpacker_outputs"
            output_dir = run_dir / "af3score_round1"
        else:
            input_pdb_dir = run_dir / "flowpacker_round2" / "flowpacker_outputs"
            output_dir = run_dir / "af3score_round2"
        cfg["command"] = [
            "python",
            str(script),
            "--input_pdb_dir",
            str(input_pdb_dir),
            "--output_dir",
            str(output_dir),
            "--af3score_repo",
            str(af3_repo),
            "--model_dir",
            str(af3_weights),
            "--num_jobs",
            str(num_jobs),
            "--write_cif_model",
            "true",
            "--export_pdb_dir",
            str(output_dir / "pdbs"),
        ]
        af3_db = tools.get("af3_db")
        if af3_db:
            cfg["command"].extend(["--db_dir", str(af3_db)])
        return

    if step_name == "af3_refold":
        af3_repo = tools.get("af3score_repo")
        af3_weights = tools.get("af3_weights")
        if not af3_repo or not af3_weights:
            return
        script = root / "scripts" / "run_af3score.py"
        if not script.exists():
            return
        run_dir = out_dir / "output"
        num_jobs = int(os.environ.get("PPIFLOW_AF3_JOBS", "1"))
        input_pdb_dir = run_dir / "relax"
        output_dir = run_dir / "af3_refold"
        refold_cfg = (input_data.get("filters") or {}).get("af3_refold") or {}
        num_samples = refold_cfg.get("num_samples")
        model_seeds = refold_cfg.get("model_seeds")
        no_templates = refold_cfg.get("no_templates", True)
        cfg["command"] = [
            "python",
            str(script),
            "--input_pdb_dir",
            str(input_pdb_dir),
            "--output_dir",
            str(output_dir),
            "--af3score_repo",
            str(af3_repo),
            "--model_dir",
            str(af3_weights),
            "--num_jobs",
            str(num_jobs),
            "--write_cif_model",
            "true",
            "--export_pdb_dir",
            str(output_dir / "pdbs"),
        ]
        if num_samples is not None:
            cfg["command"].extend(["--num_samples", str(int(num_samples))])
        if model_seeds:
            if isinstance(model_seeds, (list, tuple)):
                model_seeds = ",".join(str(s) for s in model_seeds)
            cfg["command"].extend(["--model_seeds", str(model_seeds)])
        if no_templates:
            cfg["command"].append("--no_templates")
        af3_db = tools.get("af3_db")
        if af3_db:
            cfg["command"].extend(["--db_dir", str(af3_db)])
        target_offsets = (input_data.get("target") or {}).get("chain_offsets")
        if target_offsets:
            cfg["command"].extend(["--target_offsets_json", str(target_offsets)])
        cfg["command"].extend(["--target_chain", "B"])
        return

    if step_name == "dockq":
        dockq_bin = tools.get("dockq_bin")
        if not dockq_bin:
            return
        script = root / "scripts" / "run_dockq.py"
        if not script.exists():
            return
        run_dir = out_dir / "output"
        use_refold = bool(input_data.get("optional_af3_refold") and tools.get("af3score_repo") and tools.get("af3_weights"))
        input_pdb_dir = run_dir / "af3_refold" / "pdbs" if use_refold else run_dir / "relax"
        reference_pdb_dir = run_dir / "af3score_round2" / "filtered_pdbs"
        output_dir = run_dir / "dockq"
        cfg["command"] = [
            "python",
            str(script),
            "--dockq_bin",
            str(dockq_bin),
            "--input_pdb_dir",
            str(input_pdb_dir),
            "--reference_pdb_dir",
            str(reference_pdb_dir),
            "--output_dir",
            str(output_dir),
            "--allowed_mismatches",
            "10",
            "--skip_existing",
        ]
        return

def configure_pipeline(args) -> dict:
    out_dir = Path(args.output).resolve()
    ensure_dir(out_dir)

    if args.input:
        spec = load_input(args.input)
        data = spec.data
        base_dir = Path(args.input).resolve().parent
    else:
        spec = build_input_from_cli(args)
        data = spec.data
        base_dir = Path.cwd()

    data = apply_preset(data, args.preset)
    validate_input(data)
    normalized = normalize_input(data, base_dir=base_dir, output_dir=out_dir)
    try:
        gap = int(normalized.get("target_concat_gap") or 50)
        chain_map = (normalized.get("target") or {}).get("chain_map")
        print(f"Target chains concatenated into chain B (gap={gap}).")
        if chain_map:
            print(f"Mapping written to {chain_map}.")
    except Exception:
        pass
    normalized = _auto_rename_framework(normalized, out_dir)

    # write CLI input yaml if needed
    if spec.source_path is None:
        write_cli_input_yaml(data, out_dir)

    # write pipeline_input.json (resume identity)
    input_json_path = out_dir / "pipeline_input.json"
    write_json(input_json_path, normalized, indent=2)

    tool_versions = collect_tool_versions(normalized.get("tools") or {})
    sampling = normalized.get("sampling") or {}
    target_n = int(sampling.get("samples_per_target", 0) or 0)
    state = init_or_update_state(
        out_dir=out_dir,
        input_sha256=sha256_json(normalized),
        tool_versions=tool_versions,
        target_n=target_n,
        seeds=(sampling.get("seeds") or None),
    )

    # write per-step configs
    config_dir = out_dir / "config"
    if config_dir.exists():
        counter = 1
        while config_dir.with_name(f"previous-config-{counter}").exists():
            counter += 1
        config_dir.rename(config_dir.with_name(f"previous-config-{counter}"))
    ensure_dir(config_dir)

    steps_info = []
    steps_to_write = list(STEP_ORDER)
    dockq_min = ((normalized.get("filters") or {}).get("dockq") or {}).get("min")
    if dockq_min is None:
        steps_to_write = [s for s in steps_to_write if s != "dockq"]
    if not normalized.get("optional_af3_refold"):
        steps_to_write = [s for s in steps_to_write if s != "af3_refold"]
    run_id = int((state.get("runs") or [{"run_id": 0}])[0].get("run_id", 0))
    for step_name in steps_to_write:
        cfg = _step_config(step_name, run_id, normalized)
        cfg["input"] = normalized
        _apply_default_command(step_name, run_id, out_dir, normalized, cfg)
        cfg_path = config_dir / f"step_{step_name}.yaml"
        write_yaml(cfg_path, cfg)
        steps_info.append({"name": step_name, "config_file": str(cfg_path.relative_to(out_dir))})

    steps_yaml = {"steps": steps_info}
    write_yaml(out_dir / "steps.yaml", steps_yaml)

    return {
        "out_dir": str(out_dir),
        "input": normalized,
        "state": state,
    }
