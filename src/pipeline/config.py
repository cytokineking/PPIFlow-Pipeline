from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable

from . import PIPELINE_VERSION
from .io import read_yaml, repo_root, resolve_optional_path, resolve_path, write_yaml
from .hotspots import expand_hotspots, parse_chain_list
from .target_concat import (
    concatenate_target_chains,
    compress_hotspots,
    map_hotspots_to_concatenated,
    maybe_write_hotspots_file,
)


_ROOT = repo_root()
DEFAULT_CONFIG_PATHS = {
    "binder": str(_ROOT / "src/configs/inference_binder.yaml"),
    "antibody": str(_ROOT / "src/configs/inference_nanobody.yaml"),
    "vhh": str(_ROOT / "src/configs/inference_nanobody.yaml"),
    "binder_partial": str(_ROOT / "src/configs/inference_binder_partial.yaml"),
    "antibody_partial": str(_ROOT / "src/configs/inference_nanobody.yaml"),
    "vhh_partial": str(_ROOT / "src/configs/inference_nanobody.yaml"),
}

PRESETS = {
    "fast": {
        "sampling": {"samples_per_target": 50},
        "optional_af3_refold": False,
    },
    "full": {
        "sampling": {"samples_per_target": 200},
        "optional_af3_refold": True,
    },
    "custom": {},
}


class ConfigError(ValueError):
    pass


@dataclass
class InputSpec:
    data: Dict[str, Any]
    source_path: str | None = None


def _apply_default(dst: Dict[str, Any], src: Dict[str, Any]) -> None:
    for k, v in src.items():
        if isinstance(v, dict):
            dst.setdefault(k, {})
            if isinstance(dst[k], dict):
                _apply_default(dst[k], v)
        else:
            dst.setdefault(k, v)


def apply_preset(input_data: Dict[str, Any], preset: str) -> Dict[str, Any]:
    if preset not in PRESETS:
        raise ConfigError(f"Unknown preset: {preset}")
    if preset == "custom":
        return input_data
    defaults = PRESETS[preset]
    merged = json.loads(json.dumps(input_data))
    _apply_default(merged, defaults)
    return _apply_binder_defaults(merged)


def _apply_binder_defaults(data: Dict[str, Any]) -> Dict[str, Any]:
    if data.get("protocol") != "binder":
        return data
    seq = data.setdefault("sequence_design", {})
    r1 = seq.setdefault("round1", {})
    r2 = seq.setdefault("round2", {})
    r1.setdefault("sampling_temp", 0.2)
    r1.setdefault("num_seq_per_backbone", 16)
    r1.setdefault("bias_large_residues", True)
    r1.setdefault("bias_num", 8)
    r1.setdefault("bias_residues", ["F", "M", "W"])
    r1.setdefault("bias_weight", 0.7)
    r2.setdefault("sampling_temp", 0.1)
    r2.setdefault("num_seq_per_backbone", 4)
    r2.setdefault("use_soluble_ckpt", True)
    return data


def validate_input(data: Dict[str, Any]) -> None:
    protocol = data.get("protocol")
    if protocol not in {"binder", "antibody", "vhh"}:
        raise ConfigError("protocol must be one of: binder, antibody, vhh")
    if not data.get("name"):
        raise ConfigError("name is required")
    target = data.get("target") or {}
    if not target.get("pdb"):
        raise ConfigError("target.pdb is required")
    if not target.get("chains"):
        raise ConfigError("target.chains is required")
    if protocol == "binder":
        binder = data.get("binder") or {}
        if not binder.get("length"):
            raise ConfigError("binder.length is required for binder protocol")
        if data.get("framework"):
            raise ConfigError("framework block must be absent for binder protocol")
    else:
        framework = data.get("framework") or {}
        if not framework.get("pdb"):
            raise ConfigError("framework.pdb is required for antibody/vhh")
        if not framework.get("heavy_chain"):
            raise ConfigError("framework.heavy_chain is required for antibody/vhh")
        if protocol == "antibody" and not framework.get("light_chain"):
            raise ConfigError("framework.light_chain is required for antibody")
        if protocol == "vhh" and framework.get("light_chain"):
            raise ConfigError("framework.light_chain must be omitted for vhh")
        if not framework.get("cdr_length"):
            raise ConfigError("framework.cdr_length is required for antibody/vhh")


def normalize_input(
    data: Dict[str, Any],
    *,
    base_dir: str | Path | None = None,
    output_dir: str | Path | None = None,
) -> Dict[str, Any]:
    out = _apply_binder_defaults(json.loads(json.dumps(data)))
    out["pipeline_version"] = PIPELINE_VERSION
    protocol = out.get("protocol")
    out["target"] = out.get("target") or {}
    out["target"]["pdb"] = resolve_path(out["target"]["pdb"], base_dir=base_dir)
    chains = parse_chain_list(out["target"].get("chains"))
    if chains:
        out["target"]["chains"] = chains
    out.setdefault("target_concat_gap", 50)
    out["target_concat_enabled"] = True
    if output_dir is None:
        raise ConfigError("normalize_input requires output_dir for target concatenation")
    try:
        concat_info = concatenate_target_chains(
            out["target"]["pdb"],
            out["target"]["chains"],
            output_dir,
            gap_residues=int(out.get("target_concat_gap") or 50),
        )
    except Exception as exc:
        raise ConfigError(f"Failed to concatenate target chains: {exc}") from exc
    out["target"]["original_pdb"] = out["target"]["pdb"]
    out["target"]["pdb"] = concat_info["concatenated_pdb"]
    out["target"]["chains"] = ["B"]
    out["target"]["chain_map"] = concat_info["chain_map_path"]
    out["target"]["chain_offsets"] = concat_info["offsets_path"]

    if "hotspots" in out["target"]:
        try:
            expanded = expand_hotspots(
                out["target"].get("hotspots"),
                pdb_path=out["target"]["original_pdb"],
            )
            mapped = map_hotspots_to_concatenated(expanded, concat_info["chain_map_entries"])
            compressed = compress_hotspots(mapped)
            mapped_value, mapped_file = maybe_write_hotspots_file(compressed, output_dir)
            out["target"]["hotspots"] = mapped_value
            if mapped_file:
                out["target"]["hotspots_file"] = mapped_file
        except Exception as exc:
            raise ConfigError(f"Failed to expand hotspots: {exc}") from exc
    if out.get("framework"):
        out["framework"]["pdb"] = resolve_path(out["framework"]["pdb"], base_dir=base_dir)
    tools = out.get("tools") or {}
    for key in list(tools.keys()):
        tools[key] = resolve_optional_path(tools[key], base_dir=base_dir)
    out["tools"] = tools
    if not tools.get("rosetta_bin"):
        candidate = _ROOT / "assets" / "tools" / "rosetta_scripts"
        if candidate.exists():
            tools["rosetta_bin"] = str(candidate)
    if not tools.get("rosetta_db"):
        candidate = _ROOT / "assets" / "tools" / "rosetta_db"
        if candidate.exists():
            tools["rosetta_db"] = str(candidate)
    if not tools.get("af3_db"):
        env_db = os.environ.get("AF3_DB_DIR") or os.environ.get("PPIFLOW_AF3_DB")
        candidate = None
        if env_db:
            candidate = env_db
        else:
            for path in (_ROOT / "assets" / "tools" / "af3_db", _ROOT / "assets" / "tools" / "af3_db_stub"):
                if path.exists():
                    candidate = str(path)
                    break
        if candidate:
            tools["af3_db"] = resolve_optional_path(candidate, base_dir=base_dir)
    if not tools.get("mpnn_ckpt"):
        mpnn_repo = tools.get("mpnn_repo")
        if mpnn_repo:
            candidate = Path(mpnn_repo) / "model_weights"
            if candidate.exists():
                tools["mpnn_ckpt"] = str(candidate)
    if not tools.get("abmpnn_ckpt"):
        candidate = _ROOT / "assets" / "weights" / "abmpnn" / "abmpnn.pt"
        if candidate.exists():
            tools["abmpnn_ckpt"] = str(candidate)

    # Protocol-aware defaults (paper-aligned)
    seq = out.get("sequence_design") or {}
    r1 = seq.get("round1") or {}
    r2 = seq.get("round2") or {}
    if protocol == "binder":
        r1.setdefault("sampling_temp", 0.2)
        r1.setdefault("num_seq_per_backbone", 16)
        r1.setdefault("bias_large_residues", True)
        r1.setdefault("bias_num", 8)
        r1.setdefault("bias_residues", ["F", "M", "W"])
        r1.setdefault("bias_weight", 0.7)
        r2.setdefault("sampling_temp", 0.1)
        r2.setdefault("num_seq_per_backbone", 4)
        r2.setdefault("use_soluble_ckpt", True)
    else:
        r1.setdefault("sampling_temp", 0.5)
        r1.setdefault("num_seq_per_backbone", 8)
        r1.setdefault("bias_large_residues", False)
        r1.setdefault("bias_num", 0)
        r2.setdefault("sampling_temp", 0.1)
        r2.setdefault("num_seq_per_backbone", 4)
    seq["round1"] = r1
    seq["round2"] = r2
    out["sequence_design"] = seq

    filters = out.get("filters") or {}
    af3 = filters.get("af3score") or {}
    af3.setdefault("round1", {})
    af3.setdefault("round2", {})
    af3["round1"].setdefault("iptm_min", 0.2)
    if protocol == "binder":
        af3["round1"].setdefault("ptm_min", 0.2)
    else:
        af3["round1"].setdefault("ptm_min", None)
    af3["round1"].setdefault("top_k", None)
    af3["round2"].setdefault("iptm_min", 0.5)
    af3["round2"].setdefault("ptm_min", 0.8)
    af3["round2"].setdefault("top_k", None)
    filters["af3score"] = af3
    filters.setdefault("af3_refold", {})
    filters["af3_refold"].setdefault("iptm_min", 0.7)
    filters["af3_refold"].setdefault("ptm_min", 0.8)
    filters["af3_refold"].setdefault("dockq_min", 0.49)
    filters["af3_refold"].setdefault("num_samples", 5)
    filters["af3_refold"].setdefault("model_seeds", "0-19")
    filters["af3_refold"].setdefault("no_templates", True)
    filters.setdefault("rosetta", {})
    filters["rosetta"].setdefault("interface_energy_min", -5.0)
    filters["rosetta"].setdefault("interface_distance", 12.0)
    filters["rosetta"].setdefault("relax_max_iter", 170)
    filters["rosetta"].setdefault("relax_fixbb", False)
    filters["rosetta"].setdefault("fixed_chains", "")
    filters.setdefault("dockq", {})
    filters["dockq"].setdefault("min", 0.49)
    out["filters"] = filters

    partial = out.get("partial") or {}
    partial.setdefault("start_t", 0.6)
    partial.setdefault("samples_per_target", 8)
    out["partial"] = partial

    ranking = out.get("ranking") or {}
    ranking.setdefault("top_k", 30)
    ranking.setdefault("composite_score", "iptm*100 - interface_score")
    out["ranking"] = ranking
    return out


def load_input(path: str | Path) -> InputSpec:
    data = read_yaml(path)
    if not isinstance(data, dict):
        raise ConfigError("Input YAML must be a mapping")
    return InputSpec(data=data, source_path=str(path))


def build_input_from_cli(args: Any) -> InputSpec:
    def _env_fallback(value: str | None, *keys: str) -> str | None:
        if value:
            return value
        for key in keys:
            env_value = os.environ.get(key)
            if env_value:
                return env_value
        return None

    protocol = args.protocol
    if not protocol:
        raise ConfigError("--protocol is required when no --input is provided")
    target = {
        "pdb": args.target_pdb,
        "chains": args.target_chains,
        "hotspots": args.hotspots,
    }
    if protocol == "binder":
        binder = {"length": args.binder_length}
        framework = None
    else:
        binder = None
        framework = {
            "pdb": args.framework_pdb,
            "heavy_chain": args.heavy_chain,
            "light_chain": args.light_chain,
            "cdr_length": args.cdr_length,
        }
    sampling = {"samples_per_target": args.samples_per_target}
    if protocol in {"antibody", "vhh"} and getattr(args, "vhh_backbones", None):
        sampling["samples_per_target"] = args.vhh_backbones
    sequence_design = {
        "round1": {
            "sampling_temp": args.seq1_temp,
            "num_seq_per_backbone": args.seq1_num_per_backbone
            or (getattr(args, "vhh_cdr1_num", None) if protocol in {"antibody", "vhh"} else None),
            "bias_large_residues": args.seq1_bias_large_residues if args.seq1_bias_large_residues else None,
            "bias_num": args.seq1_bias_num,
            "bias_residues": args.seq1_bias_residues,
            "bias_weight": args.seq1_bias_weight,
        },
        "round2": {
            "sampling_temp": args.seq2_temp,
            "num_seq_per_backbone": args.seq2_num_per_backbone
            or (getattr(args, "vhh_cdr2_num", None) if protocol in {"antibody", "vhh"} else None),
            "use_soluble_ckpt": args.seq2_use_soluble_ckpt if args.seq2_use_soluble_ckpt else None,
        },
    }
    filters = {
        "af3score": {
            "round1": {
                "iptm_min": args.af3score1_iptm_min,
                "ptm_min": args.af3score1_ptm_min,
                "top_k": args.af3score1_top_k,
            },
            "round2": {
                "iptm_min": args.af3score2_iptm_min,
                "ptm_min": args.af3score2_ptm_min,
                "top_k": args.af3score2_top_k,
            },
        },
        "af3_refold": {
            "iptm_min": args.af3refold_iptm_min,
            "ptm_min": args.af3refold_ptm_min,
            "dockq_min": args.af3refold_dockq_min,
            "num_samples": args.af3refold_num_samples,
            "model_seeds": args.af3refold_model_seeds,
            "no_templates": args.af3refold_no_templates,
        },
        "rosetta": {
            "interface_energy_min": args.interface_energy_min,
            "interface_distance": args.interface_distance,
            "relax_max_iter": args.relax_max_iter,
            "relax_fixbb": args.relax_fixbb,
            "fixed_chains": args.fixed_chains,
        },
        "dockq": {"min": args.dockq_min},
    }
    partial_samples = args.partial_samples_per_target
    if partial_samples is None and protocol in {"antibody", "vhh"}:
        partial_samples = getattr(args, "vhh_partial_num", None)
    partial = {
        "start_t": args.partial_start_t,
        "samples_per_target": partial_samples,
    }
    ranking = {"top_k": args.rank_top_k}
    data = {
        "protocol": protocol,
        "name": args.name,
        "target": target,
        "binder": binder,
        "framework": framework,
        "sampling": sampling,
        "sequence_design": sequence_design,
        "filters": filters,
        "partial": partial,
        "ranking": ranking,
        "optional_af3_refold": bool(args.af3_refold) if args.af3_refold else None,
        "tools": {
            "ppiflow_ckpt": args.ppiflow_ckpt,
            "abmpnn_ckpt": _env_fallback(
                args.abmpnn_ckpt,
                "ABMPNN_WEIGHTS",
                "ABMPNN_WEIGHTS_FILE",
                "PPIFLOW_ABMPNN_WEIGHTS",
            ),
            "mpnn_ckpt": _env_fallback(args.mpnn_ckpt, "MPNN_WEIGHTS", "PPIFLOW_MPNN_WEIGHTS"),
            "mpnn_ckpt_soluble": _env_fallback(
                args.mpnn_ckpt_soluble,
                "MPNN_SOLUBLE_WEIGHTS",
                "PPIFLOW_MPNN_SOLUBLE_WEIGHTS",
            ),
            "af3score_repo": _env_fallback(args.af3score_repo, "AF3SCORE_REPO", "PPIFLOW_AF3SCORE_REPO"),
            "af3_db": _env_fallback(None, "AF3_DB_DIR", "PPIFLOW_AF3_DB"),
            "rosetta_bin": args.rosetta_bin,
            "rosetta_db": args.rosetta_db,
            "flowpacker_repo": _env_fallback(args.flowpacker_repo, "FLOWPACKER_REPO", "PPIFLOW_FLOWPACKER_REPO"),
            "af3_weights": _env_fallback(args.af3_weights, "AF3_WEIGHTS", "PPIFLOW_AF3_WEIGHTS"),
            "mpnn_repo": _env_fallback(args.mpnn_repo, "PROTEINMPNN_REPO", "PPIFLOW_MPNN_REPO"),
            "abmpnn_repo": _env_fallback(args.abmpnn_repo, "PROTEINMPNN_REPO", "PPIFLOW_ABMPNN_REPO"),
            "mpnn_run": args.mpnn_run,
            "abmpnn_run": args.abmpnn_run,
            "dockq_bin": _env_fallback(args.dockq_bin, "DOCKQ_BIN", "PPIFLOW_DOCKQ_BIN"),
        },
    }
    # Strip None values
    def _strip_none(obj: Any) -> Any:
        if isinstance(obj, dict):
            return {k: _strip_none(v) for k, v in obj.items() if v is not None}
        if isinstance(obj, list):
            return [_strip_none(v) for v in obj]
        return obj

    data = _strip_none(data)
    return InputSpec(data=data, source_path=None)


def write_cli_input_yaml(input_data: Dict[str, Any], out_dir: str | Path) -> str:
    cfg_dir = Path(out_dir) / "config"
    cfg_dir.mkdir(parents=True, exist_ok=True)
    yaml_path = cfg_dir / "pipeline_input.yaml"
    write_yaml(yaml_path, input_data)
    return str(yaml_path)
