"""
Backbone Structure Sampling Pipeline - Binder Design
"""

import os
import random
import re
import dataclasses
import yaml
import time
import numpy as np
import pandas as pd
from Bio import PDB
from typing import Dict, Any, List, Optional
import argparse
import sys
from pathlib import Path
from omegaconf import OmegaConf

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from experiments.inference_antibody_partial import Experiment
from data import utils as du
from data import parsers
from data import errors
from pipeline.hotspots import expand_hotspots, parse_chain_list, resolve_hotspots_input


class ConfigManager:
    """Manages reading, modifying, and saving YAML configuration files."""

    def __init__(self, default_config_path: str):
        self.default_config = self._load_config(default_config_path)
        self.current_config: Optional[Dict[str, Any]] = None

    def _load_config(self, config_path: str) -> Dict[str, Any]:
        """Load YAML configuration file."""
        with open(config_path, "r") as f:
            return yaml.safe_load(f)

    def update_config(self, updates: Dict[str, Any]) -> Dict[str, Any]:
        """
        Deep update configuration dictionary.
        Example: updates = {'model': {'dropout': 0.2}, 'input': {'pdb_path': 'new_path.pdb'}}
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
        """Save configuration to YAML file."""
        config_to_save = self.current_config
        with open(save_path, "w") as f:
            yaml.dump(
                config_to_save,
                f,
                sort_keys=False,
                default_flow_style=False,
            )
        print(f"Config saved to: {save_path}")


# =============================================================================
# Helpers
# =============================================================================


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


def _resolve_sample_ids(sample_ids: Optional[str], samples_per_target: int) -> List[int]:
    parsed = _parse_sample_ids(sample_ids)
    if parsed is None:
        return list(range(int(samples_per_target)))
    return parsed


# =============================================================================
# Preprocessing Utilities
# =============================================================================


def extract_features(structure, struct_chains):
    """Extract features from protein structure chains.

    Args:
        structure: Biopython structure object containing the full protein complex.
        struct_chains: Dictionary mapping chain IDs to chain objects.

    Returns:
        List of dictionaries containing extracted features for each chain.
    """
    struct_feats = []
    all_seqs = set()
    for chain_id, chain in struct_chains.items():
        chain_id = du.chain_str_to_int(chain_id)
        chain_prot = parsers.process_chain(chain, chain_id)
        chain_dict = dataclasses.asdict(chain_prot)
        chain_dict = du.parse_chain_feats(
            chain_dict, normalize_positions=False
        )
        all_seqs.add(tuple(chain_dict["aatype"]))
        struct_feats.append(chain_dict)
    return struct_feats


def random_coord() -> List[float]:
    return [random.uniform(-1, 1) for _ in range(3)]


def make_virtual_ala(res_id):
    """Construct an ALA (Alanine) residue with full backbone and CB atom.

    Creates a virtual Alanine residue with random initial coordinates for all
    backbone atoms (N, CA, C, O) and the CB sidechain atom. This is used for
    initializing masked or placeholder residues.

    Args:
        res_id: The residue ID to assign to the new ALA residue.

    Returns:
        PDB.Residue.Residue object representing an Alanine with all heavy atoms.
    """
    res = PDB.Residue.Residue((" ", res_id, " "), "ALA", "")

    atom_names = ["N", "CA", "C", "O", "CB"]
    for i, atom_name in enumerate(atom_names):
        coord = random_coord()
        atom = PDB.Atom.Atom(
            atom_name,
            coord,
            1.0,
            1.0,
            " ",
            atom_name,
            i,
            element=atom_name[0],
        )
        res.add(atom)
    return res


# =============================================================================
# Index Parsing Utilities
# =============================================================================


def get_indices_from_spec(structure, spec) -> List[int]:
    """Parse global residue indices from expanded hotspot list."""
    if not spec:
        return []
    if isinstance(spec, str):
        target_ids = set([t.strip() for t in spec.split(",") if t.strip()])
    else:
        target_ids = set(str(t).strip() for t in spec if str(t).strip())
    indices = []
    res_idx = 0
    for chain in structure.get_chains():
        cid = chain.id
        for res in chain:
            if PDB.is_aa(res, standard=True):
                res_num = str(res.get_id()[1])
                icode = res.get_id()[2].strip()
                full_id = f"{cid}{res_num}{icode}"
                if full_id in target_ids:
                    indices.append(res_idx)
                res_idx += 1
    return indices


def expand_ranges(s: str | None) -> str:
    """Expand numeric range specifications in a residue string.

    Converts compact range notation like "H50-55" into individual residue
    specifications "H50,H51,H52,H53,H54,H55". Mixed formats (ranges and
    individual residues) are handled correctly.

    Args:
        s: String containing residue specifications, possibly with ranges.
            Example: "H50-55,H60,H62-65"

    Returns:
        String with all ranges expanded to individual residue specifications.
        Returns original string if no ranges are present.
    """
    if not s:
        return ""
    if "-" not in s:
        return s

    result = []
    parts = s.split(",")

    for part in parts:
        part = part.strip()
        match = re.match(r"([A-Za-z]+)(\d+)-(\d+)", part)

        if match:
            prefix, start, end = match.groups()
            start, end = int(start), int(end)
            expanded = [f"{prefix}{i}" for i in range(start, end + 1)]
            result.extend(expanded)
        else:
            result.append(part)

    return ",".join(result)


def process_file(input_info, write_dir, id):
    pdb_name = input_info["pdb_name"]

    filepath = input_info["complex_pdb"]
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(f"{pdb_name}_sample_{id}", filepath)
    res_counts = {}
    for chain in structure.get_chains():
        chain_id = chain.id
        count = sum(1 for r in chain if PDB.is_aa(r, standard=True))
        res_counts[chain_id] = count

    if input_info["light_chain"] is None:
        ab_chains = [input_info["heavy_chain"]]
    else:
        ab_chains = [input_info["heavy_chain"], input_info["light_chain"]]

    ag_chains = list(input_info["antigen_chain"])

    ab_length = sum(res_counts[ab_chain] for ab_chain in ab_chains)
    ag_length = sum(res_counts[ag_chain] for ag_chain in ag_chains)
    total_length = ab_length + ag_length

    cdr_position = input_info.get("cdr_position")
    if cdr_position:
        cdr_position_trans = expand_ranges(cdr_position)
    else:
        # Default to entire binder chain(s) when CDR positions are not provided.
        cdr_position_trans = expand_hotspots(ab_chains, pdb_path=filepath)
    cdr_indices = get_indices_from_spec(structure, cdr_position_trans)
    cdr_mask = np.zeros(total_length, dtype=np.int8)
    cdr_mask[cdr_indices] = 1

    fix_structure = input_info["fixed_positions"]
    fix_structure_trans = expand_ranges(fix_structure)
    fix_structure_indices = get_indices_from_spec(
        structure, fix_structure_trans
    )
    fix_structure_mask = np.zeros(total_length, dtype=np.int8)
    fix_structure_mask[ab_length:] = 1
    fix_structure_mask[fix_structure_indices] = 1

    hotspot_spec = expand_hotspots(
        input_info.get("hotspots"),
        pdb_path=filepath,
    )
    hotspot_index = get_indices_from_spec(structure, hotspot_spec)

    metadata = {}
    metadata["id"] = f"{pdb_name}_sample_{id}"
    metadata["pdb_name"] = f"{pdb_name}_sample_{id}"
    metadata["sample_id"] = int(id)
    processed_path = os.path.join(write_dir, f"{pdb_name}_sample_{id}.pkl")

    metadata.update(
        {
            "target_id": input_info["antigen_chain"],
            "binder_id": ab_chains,
            "processed_path": processed_path,
            "data_level": 1,
            "hotspots": hotspot_index,
        }
    )

    struct_chains = {
        chain.id.upper(): chain for chain in structure.get_chains()
    }
    struct_feats = extract_features(structure, struct_chains)
    complex_feats = du.concat_np_features(struct_feats, False)

    metadata["num_chains"] = len(struct_chains)

    # Process geometry features
    complex_aatype = complex_feats["aatype"]
    metadata["seq_len"] = len(complex_aatype)
    modeled_idx = np.where(complex_aatype != 20)[0]
    if np.sum(complex_aatype != 20) == 0:
        raise errors.LengthError("No modeled residues")
    metadata["modeled_seq_len"] = len(modeled_idx)
    complex_feats["modeled_idx"] = modeled_idx

    # chain_groups (0=Target/Antigen, 1=Binder/Antibody)
    ag_chain_ints = [
        du.chain_str_to_int(c) for c in input_info["antigen_chain"]
    ]
    ab_chain_ints = [du.chain_str_to_int(c) for c in ab_chains]

    complex_feats["chain_groups"] = np.where(
        np.isin(complex_feats["chain_index"], ag_chain_ints),
        0,
        np.where(
            np.isin(complex_feats["chain_index"], ab_chain_ints),
            1,
            complex_feats["chain_index"],
        ),
    )

    complex_feats["fix_structure_mask"] = fix_structure_mask
    complex_feats["cdr_mask"] = cdr_mask

    du.write_pkl(processed_path, complex_feats)

    return metadata


def preprocess_csv_and_pkl(args, output_dir) -> str:
    csv_path = os.path.join(output_dir, f"{args.name}_input.csv")

    antigen_chains = parse_chain_list(args.antigen_chain)
    hotspots_spec = resolve_hotspots_input(
        args.specified_hotspots, hotspots_file=args.hotspots_file
    )
    input_info = {
        "antigen_chain": antigen_chains,
        "heavy_chain": args.heavy_chain,
        "light_chain": args.light_chain,
        "pdb_name": args.name,
        "hotspots": hotspots_spec,
        "fixed_positions": args.fixed_positions,
        "cdr_position": args.cdr_position,
        "complex_pdb": args.complex_pdb,  # Pass new complex_pdb arg
    }

    sample_ids = _resolve_sample_ids(args.sample_ids, args.samples_per_target)
    for idx, id in enumerate(sample_ids):
        metadata = process_file(input_info, write_dir=output_dir, id=id)
        metadata_df = pd.DataFrame([metadata])
        header = False if idx > 0 else True
        metadata_df.to_csv(csv_path, index=False, mode="a", header=header)

    return csv_path


# =============================================================================
# Main Pipeline
# =============================================================================


def run_pipeline(args):
    """Execute the full antibody backbone sampling pipeline.

    This is the main entry point for the binder design pipeline. It orchestrates:
    1. Data preprocessing (PDB parsing, feature extraction, mask creation)
    2. Configuration management (YAML config loading and updates)
    3. Model inference (running the diffusion model for structure generation)

    Args:
        args: Command-line arguments namespace containing all pipeline parameters.
    """
    print(f"Running pipeline with args: {args}")
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)

    # Step 1: Preprocessing - Parse PDB and extract features
    input_data_dir = os.path.join(output_dir, "input")
    os.makedirs(input_data_dir, exist_ok=True)
    processed_csv_path = preprocess_csv_and_pkl(
        args=args, output_dir=input_data_dir
    )

    # Step 2: Configuration management - Load and update YAML config
    conf = ConfigManager(args.config)
    ppi_dataset_cfg: Dict[str, Any] = {
        "test_csv_path": processed_csv_path,
        "samples_per_target": 1,
        "define_hotspots": True,
    }
    if args.seed is not None:
        ppi_dataset_cfg["seed"] = int(args.seed)
    sample_ids_override = _parse_sample_ids(args.sample_ids)
    if sample_ids_override is not None:
        ppi_dataset_cfg["sample_ids_override"] = sample_ids_override

    update_configs = {
        "ppi_dataset": ppi_dataset_cfg,
        "experiment": {
            "testing_model": {
                "ckpt_path": args.model_weights,
                "save_dir": args.output_dir,
            },
            "checkpointer": {"dirpath": args.output_dir},
            "retry_Limit": args.retry_Limit,
        },
        "interpolant": {
            "min_t": args.start_t if args.start_t else 0.8,
        },
    }
    conf.update_config(update_configs)

    os.makedirs(f"{args.output_dir}/yaml", exist_ok=True)
    config_save_path = os.path.join(
        output_dir, f"yaml/antibody_sample_config_{args.name}.yml"
    )
    conf.save_config(config_save_path)

    with open(config_save_path, "r") as f:
        cfg_dict = yaml.safe_load(f)
    cfg = OmegaConf.create(cfg_dict)

    print("\nInitializing model...")
    t1 = time.time()
    exp = Experiment(cfg=cfg)

    print("\nRunning inference...")
    exp.test()

    t2 = time.time()
    print(f"sample finished in {round(t2 - t1, 2)} seconds")


# =============================================================================
# Command-Line Interface
# =============================================================================


def get_parser():
    """Create argument parser for command-line interface.

    Returns:
        argparse.ArgumentParser configured with all pipeline arguments.
    """
    parser = argparse.ArgumentParser(
        description="Protein Structure Prediction Pipeline"
    )

    # Mode 2: Complex Input (New)
    parser.add_argument(
        "--complex_pdb",
        type=str,
        help="Input complex PDB file path (supersedes antigen_pdb/framework_pdb)",
    )
    parser.add_argument(
        "--fixed_positions",
        type=str,
        required=True,
        help="Key residues to fix in complex_pdb. Format: 'H26,H27,H28,L50-63' (chain ID + residue number, '-' for ranges).",
    )
    parser.add_argument(
        "--cdr_position",
        type=str,
        help="Specify CDR residues, e.g., 'H26-32,H45-56,H97-113'",
    )

    # Common Inputs
    parser.add_argument(
        "--antigen_chain",
        "--target_chain",
        "--target_chains",
        dest="antigen_chain",
        type=str,
        required=True,
        help="Chain id(s) of the antigen/target (comma-separated or concatenated)",
    )
    parser.add_argument(
        "--heavy_chain",
        type=str,
        required=True,
        help="Chain id of the heavy chain",
    )
    parser.add_argument(
        "--light_chain",
        type=str,
        default=None,
        help="Chain id of the light chain (optional, for antibody)",
    )

    parser.add_argument(
        "--specified_hotspots",
        type=str,
        help="Hotspot spec, e.g., 'C' (whole chain) or 'C62-73,C80'",
    )
    parser.add_argument(
        "--hotspots_file",
        type=str,
        help="Path to hotspots file (one per line or comma-separated).",
    )

    # Sampling Control
    parser.add_argument(
        "--start_t",
        type=float,
        required=True,
        help="Starting t value for sampling (0.0 - 1.0).",
    )
    parser.add_argument(
        "--samples_per_target",
        type=int,
        required=True,
        help="Number of generated samples per target (e.g., 100).",
    )
    parser.add_argument(
        "--sample_ids",
        type=str,
        help="Comma-separated sample ids to generate (overrides samples_per_target).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        help="Random seed for reproducibility.",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Output directory to save results.",
    )
    parser.add_argument(
        "--retry_Limit",
        type=int,
        default=10,
        help="Maximum retry attempts if sampling fails (default: 10).",
    )

    # System
    parser.add_argument(
        "--config",
        type=str,
        default=str(SRC / "configs/inference_nanobody.yaml"),
        help="Path to configuration YAML file.",
    )
    parser.add_argument(
        "--model_weights",
        type=str,
        required=True,
        help="Path to pretrained model weights file.",
    )
    parser.add_argument(
        "--name",
        type=str,
        default="flow",
        help="Name of the test target (default: 'flow').",
    )

    return parser


if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()

    run_pipeline(args)
