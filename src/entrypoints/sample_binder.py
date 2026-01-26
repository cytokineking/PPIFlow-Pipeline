"""
Backbone Structure Sampling Pipeline - Binder Design
"""

import os
import shutil
import argparse
import json
import re
import random
import sys
from pathlib import Path
from typing import Dict, Any, Optional, List
import pandas as pd
import yaml
import numpy as np
from omegaconf import OmegaConf

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from experiments.inference_binder import Experiment
from preprocessing.get_interface_util import (
    get_residue_pairs_within_distance,
)
from preprocessing.process_pdb_for_inputs import process_file
from pipeline.hotspots import expand_hotspots, parse_chain_list, resolve_hotspots_input


class ConfigManager:
    """Manages YAML configuration loading, updating, and saving."""

    def __init__(self, default_config_path: str):
        self.default_config = self._load_config(default_config_path)
        self.current_config: Optional[Dict[str, Any]] = None

    def _load_config(self, config_path: str) -> Dict[str, Any]:
        """Load YAML configuration file."""
        with open(config_path, "r") as f:
            return yaml.safe_load(f)

    def update_config(self, updates: Dict[str, Any]) -> Dict[str, Any]:
        """Deep update configuration dictionary.

        Args:
            updates: Nested dictionary with updates to apply.
                     Example: {'model': {'dropout': 0.2}}

        Returns:
            Updated configuration dictionary.
        """

        def deep_update(original: Dict, update: Dict) -> Dict:
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

    def save_config(self, save_path: str) -> None:
        """Save current configuration to YAML file."""
        with open(save_path, "w") as f:
            yaml.dump(
                self.current_config,
                f,
                sort_keys=False,
                default_flow_style=False,
            )
        print(f"Config saved to: {save_path}")


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


def preprocess_csv_and_pkl(pdb_path: str, output_dir: str, args) -> str:
    """Process PDB file to generate pkl and metadata.

    Args:
        pdb_path: Path to the input PDB file.
        output_dir: Directory to write processed files.
        args: Command-line arguments containing chain info and hotspots.

    Returns:
        Path to the generated CSV file.
    """
    print(f"Preprocessing data: {pdb_path}")
    input_info = {"pdbfile": pdb_path, "PDBID": args.name}

    target_chains = parse_chain_list(args.target_chain)
    if len(target_chains) != 1:
        raise ValueError("Binder protocol currently supports a single target chain")
    target_chain = target_chains[0]

    hotspots_spec = resolve_hotspots_input(
        args.specified_hotspots, hotspots_file=args.hotspots_file
    )
    if hotspots_spec is not None:
        print(f"Using specified hotspots: {hotspots_spec}")
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
            "Hotspots not specified, generating hotspots from binder interface..."
        )
        assert (
            args.binder_chain is not None
        ), "Binder chain required when no hotspots given"
        interface_data = get_residue_pairs_within_distance(
            pdb_path,
            target_chain,
            args.binder_chain,
            distance_threshold=10.0,
        )
        input_info["contact_pairs10A"] = interface_data[0]
        input_info["chain1_residues"] = interface_data[1]
        input_info["chain2_residues"] = interface_data[2]
        input_info["chain1_id"] = target_chain
        input_info["chain2_id"] = args.binder_chain

    metadata = process_file(input_info, write_dir=output_dir)
    metadata["num_chains"] = 2
    metadata_df = pd.DataFrame([metadata])
    csv_path = os.path.join(output_dir, f"{args.name}_input.csv")
    metadata_df.to_csv(csv_path, index=False)

    return csv_path


def run_pipeline(args) -> None:
    """Execute the complete sampling pipeline.

    Args:
        args: Parsed command-line arguments.
    """
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)

    sample_ids_override = _parse_sample_ids(args.sample_ids)
    length_schedule = _load_length_schedule(
        args.length_schedule_path, args.samples_per_target
    )

    # Preprocessing: generate CSV from PDB or use provided CSV
    if args.input_pdb is not None:
        assert args.target_chain is not None, "Target chain ID required"
        input_data_dir = os.path.join(output_dir, "input")
        os.makedirs(input_data_dir, exist_ok=True)
        processed_csv_path = preprocess_csv_and_pkl(
            pdb_path=args.input_pdb, output_dir=input_data_dir, args=args
        )
    else:
        assert (
            args.input_csv is not None
        ), "Either input_csv or input_pdb required"
        processed_csv_path = args.input_csv
        shutil.copy(
            processed_csv_path,
            os.path.join(output_dir, f"{args.name}_input.csv"),
        )

    # Load and update configuration
    config_manager = ConfigManager(args.config)
    hotspots_spec = resolve_hotspots_input(
        args.specified_hotspots, hotspots_file=args.hotspots_file
    )
    has_specified_hotspots = hotspots_spec is not None
    ppi_dataset_cfg: Dict[str, Any] = {
        "test_csv_path": processed_csv_path,
        "samples_min_length": args.samples_min_length,
        "samples_max_length": args.samples_max_length,
        "samples_per_target": args.samples_per_target,
        "define_hotspots": has_specified_hotspots,
        "min_hotspot_ratio": args.sample_hotspot_rate_min,
        "max_hotspot_ratio": args.sample_hotspot_rate_max,
    }
    if args.seed is not None:
        ppi_dataset_cfg["seed"] = int(args.seed)
    if sample_ids_override is not None:
        ppi_dataset_cfg["sample_ids_override"] = sample_ids_override
    if length_schedule is not None:
        ppi_dataset_cfg["length_schedule"] = length_schedule

    update_configs = {
        "ppi_dataset": ppi_dataset_cfg,
        "experiment": {
            "testing_model": {
                "ckpt_path": args.model_weights,
                "save_dir": output_dir,
            },
            "checkpointer": {"dirpath": output_dir},
        },
    }
    config_manager.update_config(update_configs)

    # Save config and load as OmegaConf
    config_save_path = os.path.join(output_dir, "sample_config.yml")
    config_manager.save_config(config_save_path)
    cfg = OmegaConf.create(config_manager.current_config)

    # Initialize model and run inference
    print("\nInitializing model...")
    exp = Experiment(cfg=cfg)

    print("\nRunning inference...")
    exp.test()
    print("Sampling finished.")


def get_parser() -> argparse.ArgumentParser:
    """Create and return the argument parser."""
    parser = argparse.ArgumentParser(
        description="Protein Structure Prediction Pipeline - Binder Design"
    )

    # Input options (mutually exclusive, one required)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--input_pdb", type=str, help="Input protein PDB file path"
    )
    input_group.add_argument(
        "--input_csv", type=str, help="Input CSV file path"
    )

    # Chain configuration
    parser.add_argument(
        "--target_chain",
        type=str,
        default="R",
        help="Chain ID of target protein (single chain only)",
    )
    parser.add_argument(
        "--binder_chain",
        type=str,
        default=None,
        help="Chain ID of binder protein",
    )

    # Model configuration
    parser.add_argument(
        "--config",
        type=str,
        default=str(SRC / "configs/inference_binder.yaml"),
        help="Default configuration file path",
    )

    # Hotspot configuration
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
        help="Minimum hotspot sampling rate",
    )
    parser.add_argument(
        "--sample_hotspot_rate_max",
        type=float,
        default=0.5,
        help="Maximum hotspot sampling rate",
    )

    # Sampling configuration
    parser.add_argument(
        "--samples_min_length",
        type=int,
        default=50,
        help="Minimum number of residues per sample",
    )
    parser.add_argument(
        "--samples_max_length",
        type=int,
        default=100,
        help="Maximum number of residues per sample",
    )
    parser.add_argument(
        "--samples_per_target",
        type=int,
        default=100,
        help="Number of samples to generate per target",
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

    # Checkpoint and output
    parser.add_argument(
        "--model_weights",
        type=str,
        required=True,
        help="Model weights checkpoint file path",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Output directory for results",
    )
    parser.add_argument(
        "--name",
        type=str,
        default="test_target",
        help="Name identifier for this run",
    )

    return parser


def validate_inputs(args) -> None:
    """Validate that required input files exist."""
    assert os.path.exists(
        args.config
    ), f"Config file not found: {args.config}"
    assert os.path.exists(
        args.model_weights
    ), f"Model weights not found: {args.model_weights}"
    if args.length_schedule_path is not None:
        assert os.path.exists(
            args.length_schedule_path
        ), f"Length schedule file not found: {args.length_schedule_path}"

    if args.input_pdb is not None:
        assert os.path.exists(
            args.input_pdb
        ), f"PDB file not found: {args.input_pdb}"

    target_chains = parse_chain_list(args.target_chain)
    assert len(target_chains) == 1, "Binder protocol currently supports a single target chain"
    target_chain = target_chains[0]
    hotspots_spec = resolve_hotspots_input(
        args.specified_hotspots, hotspots_file=args.hotspots_file
    )
    if hotspots_spec is not None and args.input_pdb is not None:
        expanded = expand_hotspots(hotspots_spec, pdb_path=args.input_pdb)
        bad = [h for h in expanded if h and h[0] != target_chain]
        assert not bad, "Hotspot chain ID must match target chain ID"

def main() -> None:
    """Main entry point for the binder sampling pipeline."""
    parser = get_parser()
    args = parser.parse_args()

    validate_inputs(args)

    # Hotspot rate limits only apply when hotspots are NOT specified
    hotspots_spec = resolve_hotspots_input(
        args.specified_hotspots, hotspots_file=args.hotspots_file
    )
    if hotspots_spec is None:
        if (
            args.sample_hotspot_rate_min != 0.2
            or args.sample_hotspot_rate_max != 0.5
        ):
            parser.error(
                "--specified_hotspots or --hotspots_file required when customizing sampling rates"
            )

    run_pipeline(args)


if __name__ == "__main__":
    main()
