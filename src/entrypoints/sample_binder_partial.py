"""
Backbone Structure Sampling Pipeline - Binder Design
"""

import os
import shutil
import pandas as pd
import time
import yaml
import json
import re
import random
import sys
from pathlib import Path
from typing import Dict, Any, Optional, List
import argparse
import numpy as np
from omegaconf import OmegaConf

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from experiments.inference_binder_partial import Experiment
from preprocessing.get_interface_util import (
    get_residue_pairs_within_distance,
)
from preprocessing.process_partial import process_file
from pipeline.hotspots import expand_hotspots, parse_chain_list, resolve_hotspots_input

"""
Config Manager: Read, modify and save YAML configuration files
"""


class ConfigManager:
    def __init__(self, default_config_path: str):
        self.default_config = self._load_config(default_config_path)
        self.current_config = None

    def _load_config(self, config_path: str) -> Dict[str, Any]:
        """Load YAML configuration file"""
        with open(config_path, "r") as f:
            return yaml.safe_load(f)

    def update_config(self, updates: Dict[str, Any]) -> Dict[str, Any]:
        """
        Deeply update configuration dictionary
        Example: updates = {'model': {'dropout': 0.2}, 'input': {'pdb_path': 'new_path.pdb'}}
        """

        def deep_update(original, update):
            for key, value in update.items():
                if isinstance(value, dict) and key in original:
                    deep_update(original[key], value)
                else:
                    original[key] = value
            return original

        self.current_config = deep_update(
            self.default_config.copy(), updates
        )
        return self.current_config

    def save_config(self, save_path: str):
        """Save configuration to YAML file"""
        config_to_save = self.current_config
        with open(save_path, "w") as f:
            yaml.dump(
                config_to_save,
                f,
                sort_keys=False,
                default_flow_style=False,
            )
        print(f"Config saved to: {save_path}")


# --------------------------
# Helpers
# --------------------------
def _parse_sample_ids(sample_ids: Optional[str]) -> Optional[List[int]]:
    if not sample_ids:
        return None
    ids: List[int] = []
    for token in sample_ids.split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            start, end = token.split("-", 1)
            start_i = int(start)
            end_i = int(end)
            if end_i < start_i:
                raise ValueError("sample_ids range must be ascending")
            ids.extend(list(range(start_i, end_i + 1)))
        else:
            ids.append(int(token))
    if not ids:
        raise ValueError("sample_ids is empty")
    return ids


def _load_length_schedule(path: Optional[str], expected_len: int) -> Optional[List[int]]:
    if not path:
        return None
    with open(path, "r") as f:
        data = json.load(f)
    lengths = data.get("lengths") if isinstance(data, dict) else data
    if not isinstance(lengths, list):
        raise ValueError("length_schedule_path must be a JSON list or {\"lengths\": [...]} object")
    schedule = [int(x) for x in lengths]
    if expected_len and len(schedule) != int(expected_len):
        raise ValueError("length_schedule length must equal samples_per_target")
    return schedule


# --------------------------
# Preprocessing Module
# --------------------------
def get_motif_residues(motif_contig, chainid):
    """
    param:
    motif_contig (str): String containing ranges, e.g., "L1-3,L10,L12-13".

    return:
    list: Nested list, e.g., [[1, 2, 3], [10], [12, 13]].
    """
    motif_contig = motif_contig.replace(chainid, "")
    ranges = motif_contig.split(",")

    result = []
    for r in ranges:
        if "-" in r:
            # If range, e.g., '1-3'
            start, end = map(int, r.split("-"))
            result.append(list(range(start, end + 1)))
        else:
            # If single number, e.g., '10'
            result.append([int(r)])
    return result


def _normalize_motif_contig(motif_contig):
    if motif_contig is None:
        return None
    if isinstance(motif_contig, float) and np.isnan(motif_contig):
        return None
    if isinstance(motif_contig, str):
        cleaned = motif_contig.strip()
        if not cleaned or cleaned.lower() in {"nan", "none", "null"}:
            return None
        return cleaned
    return motif_contig


def preprocess_csv_and_pkl(pdb_path, output_dir, args) -> str:
    """
    Process PDB file to generate pkl and metadata
    Returns: csv_path
    """
    print("preprocessing data: ", pdb_path)
    input_info = {"pdbfile": pdb_path, "PDBID": args.name}
    target_chains = parse_chain_list(args.target_chain)
    if len(target_chains) != 1:
        raise ValueError("Binder protocol currently supports a single target chain")
    target_chain = target_chains[0]
    hotspots_spec = resolve_hotspots_input(
        args.specified_hotspots, hotspots_file=args.hotspots_file
    )
    if hotspots_spec is not None:
        print("use specified hotspots: ", hotspots_spec)
        expanded = expand_hotspots(hotspots_spec, pdb_path=pdb_path)
        hotspot_residues = []
        for token in expanded:
            if not token or token[0] != target_chain:
                raise ValueError("Hotspots must be on the target chain")
            match = re.match(r"^%s(\d+)" % re.escape(target_chain), token)
            if not match:
                raise ValueError(f"Invalid hotspot token: {token}")
            hotspot_residues.append(int(match.group(1)))
        input_info["chain1_residues"] = hotspot_residues
        input_info["chain1_id"] = target_chain
        if args.binder_chain is not None:
            input_info["chain2_id"] = args.binder_chain
    else:
        print(
            "hotspots not specified, generate hotspots according to known binder interface..."
        )
        assert (
            args.binder_chain is not None
        ), "must provide binder chain if no hotspots are given"
        interface_data = get_residue_pairs_within_distance(
            pdb_path,
            target_chain,
            args.binder_chain,
            distance_threshold=args.interface_dist,
        )

        input_info["contact_pairs10A"] = interface_data[0]
        input_info["chain1_residues"] = interface_data[1]
        input_info["chain2_residues"] = interface_data[2]
        input_info["chain1_id"] = target_chain
        input_info["chain2_id"] = args.binder_chain

    motif_contig = _normalize_motif_contig(args.motif_contig)
    if motif_contig is not None:
        input_info["binder_motif"] = get_motif_residues(
            motif_contig, args.binder_chain
        )

    metadata = process_file(input_info, write_dir=output_dir)
    metadata["num_chains"] = 2
    # metadata['contig'] = args.motif_contig
    # metadata['sample_binder_len'] = args.sample_binder_len

    metadata_df = pd.DataFrame([metadata])  # one item only
    csv_path = os.path.join(output_dir, f"{args.name}_input.csv")
    metadata_df.to_csv(csv_path, index=False)

    return csv_path


# --------------------------
# Main Pipeline
# --------------------------
def run_pipeline(args):
    # 1. Create output directory
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)

    sample_ids_override = _parse_sample_ids(args.sample_ids)
    length_schedule = _load_length_schedule(
        args.length_schedule_path, args.samples_per_target
    )

    # 2. Preprocessing
    if args.input_pdb is not None:
        # print(f"\nPreprocessing {args.input_pdb} for metadata...")
        assert (
            args.target_chain is not None
        ), "must provide target chain id"
        input_data_dir = os.path.join(output_dir, "input")
        os.makedirs(input_data_dir, exist_ok=True)
        processed_csv_path = preprocess_csv_and_pkl(
            pdb_path=args.input_pdb, output_dir=input_data_dir, args=args
        )
    else:
        assert (
            args.input_csv is not None
        ), "must provide either input csv or pdb"
        processed_csv_path = args.input_csv
        shutil.copy(
            processed_csv_path,
            os.path.join(output_dir, f"{args.name}_input.csv"),
        )

    # 3. Load and update config
    conf = ConfigManager(args.config)  # load default config
    hotspots_spec = resolve_hotspots_input(
        args.specified_hotspots, hotspots_file=args.hotspots_file
    )
    motif_contig = _normalize_motif_contig(args.motif_contig)
    ppi_dataset_cfg: Dict[str, Any] = {
        "test_csv_path": processed_csv_path,
        "samples_per_target": args.samples_per_target,
        "define_hotspots": (True if (hotspots_spec is not None) else False),
        "min_hotspot_ratio": args.sample_hotspot_rate_min,
        "max_hotspot_ratio": args.sample_hotspot_rate_max,
        "motif": (
            None
            if (motif_contig is None)
            else {"define_motif": True}
        ),
    }
    if args.seed is not None:
        ppi_dataset_cfg["seed"] = int(args.seed)
    if sample_ids_override is not None:
        ppi_dataset_cfg["sample_ids_override"] = sample_ids_override
    if length_schedule is not None:
        ppi_dataset_cfg["length_schedule"] = length_schedule

    update_configs = {
        "data": {"task": "binder_motif_partial"},
        "ppi_dataset": ppi_dataset_cfg,
        "experiment": {
            "testing_model": {
                "ckpt_path": args.model_weights,
                "save_dir": args.output_dir,
            },
            "checkpointer": {"dirpath": args.output_dir},
        },
        "interpolant": {
            "min_t": args.start_t,
            "max_t": args.start_t + 0.02,
        },
    }
    conf.update_config(update_configs)
    os.makedirs(f"{args.output_dir}/yaml", exist_ok=True)
    config_save_path = os.path.join(
        output_dir, f"yaml/binder_partial_sample_config_{args.name}.yml"
    )
    conf.save_config(config_save_path)

    with open(config_save_path, "r") as f:
        cfg_dict = yaml.safe_load(f)
    cfg = OmegaConf.create(cfg_dict)

    # 4. Initialize model
    print("\nInitializing model...")
    t1 = time.time()
    exp = Experiment(cfg=cfg)

    # 5. Run inference
    print("\nRunning inference...")
    exp.test()

    t2 = time.time()
    print(f"sample finished in {round(t2 - t1, 2)} seconds")

    return


# --------------------------
# Command Line Interface
# --------------------------
def get_parser():
    parser = argparse.ArgumentParser(
        description="Protein Structure Prediction Pipeline"
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--input_pdb", type=str, help="Input protein PDB file path"
    )
    group.add_argument("--input_csv", type=str, help="Input CSV file path")

    parser.add_argument(
        "--target_chain",
        type=str,
        default="R",
        help="Chain id of the target (only single chain target is supported), e.g., B",
    )
    parser.add_argument(
        "--binder_chain",
        type=str,
        default="L",
        help="Chain id of the binder (motif), e.g., L",
    )

    # Model configs
    parser.add_argument(
        "--config",
        type=str,
        default=str(SRC / "configs/inference_binder_partial.yaml"),
        help="default configuration file",
    )

    # Hotspots
    parser.add_argument(
        "--specified_hotspots",
        type=str,
        help="Hotspot spec, e.g., 'R' (whole chain) or 'R62-73,R80'",
    )
    parser.add_argument(
        "--hotspots_file",
        type=str,
        help="Path to hotspots file (one per line or comma-separated).",
    )
    parser.add_argument(
        "--sample_hotspot_rate_min",
        type=float,
        default=0.2,
        help="Minimum sampling rate (default: 0.2)",
    )
    parser.add_argument(
        "--sample_hotspot_rate_max",
        type=float,
        default=0.5,
        help="Maximum sampling rate (default: 0.5)",
    )
    parser.add_argument(
        "--interface_dist",
        type=float,
        default=6.0,
        help="interface distance between target and binder",
    )

    # Motif
    parser.add_argument(
        "--motif_contig", type=str, help="Motif contig, e.g., 'L19-27,L31'"
    )
    
    # Samples
    parser.add_argument(
        "--samples_per_target",
        type=int,
        default=100,
        help="number of samples",
    )
    parser.add_argument(
        "--sample_ids",
        type=str,
        help="Comma-separated sample ids to generate (overrides samples_per_target)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        help="Random seed for reproducibility",
    )
    parser.add_argument(
        "--length_schedule_path",
        type=str,
        help="Path to JSON length schedule (list or {\"lengths\": [...]})",
    )

    # Checkpoint
    parser.add_argument(
        "--model_weights",
        type=str,
        default="binder.ckpt",
        help="Model weights file path",
    )

    # Interpolant t  #partial
    parser.add_argument(
        "--start_t", type=float, default=0.15, help="start_t"
    )

    # Output dir
    parser.add_argument(
        "--output_dir", type=str, required=True, help="Output directory"
    )
    parser.add_argument(
        "--name", type=str, default="test_target", help="Test target name"
    )
    return parser


if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()

    # Validate input files
    assert os.path.exists(
        args.input_pdb
    ), f"PDB file not found: {args.input_pdb}"
    assert os.path.exists(
        args.config
    ), f"Config file not found: {args.config}"
    assert os.path.exists(
        args.model_weights
    ), f"Model weights file not found: {args.model_weights}"
    if args.length_schedule_path is not None:
        assert os.path.exists(
            args.length_schedule_path
        ), f"Length schedule file not found: {args.length_schedule_path}"

    target_chains = parse_chain_list(args.target_chain)
    assert len(target_chains) == 1, "Binder protocol currently supports a single target chain"
    target_chain = target_chains[0]

    hotspots_spec = resolve_hotspots_input(
        args.specified_hotspots, hotspots_file=args.hotspots_file
    )
    if hotspots_spec:
        if (
            args.sample_hotspot_rate_min != 0.2
            or args.sample_hotspot_rate_max != 0.5
        ):
            parser.error(
                "--specified_hotspots/--hotspots_file cannot be used together with sampling rate args"
            )
        if args.input_pdb:
            expanded = expand_hotspots(hotspots_spec, pdb_path=args.input_pdb)
            bad = [h for h in expanded if h and h[0] != target_chain]
            if bad:
                parser.error("Hotspots must be on the target chain")

    os.makedirs(args.output_dir, exist_ok=True)

    run_pipeline(args)
