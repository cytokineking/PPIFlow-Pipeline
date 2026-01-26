from __future__ import annotations

import argparse
from pathlib import Path

from .configure import configure_pipeline
from .execute import execute_pipeline


def _add_common_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--protocol", type=str, choices=["binder", "antibody", "vhh"], help="Protocol")
    parser.add_argument("--preset", type=str, choices=["fast", "full", "custom"], default="full")
    parser.add_argument("--input", type=str, help="Path to design.yaml")
    parser.add_argument("--output", type=str, required=True, help="Output directory")
    parser.add_argument(
        "--steps",
        type=str,
        default="all",
        help="Step groups or step names (all|gen,seq,score,rosetta,partial,rank)",
    )
    parser.add_argument("--reuse", action="store_true", help="Reuse existing outputs")
    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Continue to next step if a step fails",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Stream step subprocess output to console (also logs to file)",
    )
    parser.add_argument("--num-workers", type=int, default=8)
    parser.add_argument("--devices", type=str, default=None)


def _add_cli_input_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--name", type=str, help="Design name")
    parser.add_argument("--target_pdb", type=str, help="Target PDB path")
    parser.add_argument("--target_chains", type=str, help="Target chain(s) (comma-separated)")
    parser.add_argument("--hotspots", type=str, help="Hotspots (comma-separated)")
    parser.add_argument("--binder_length", type=str, help="Binder length (e.g. 75-90)")
    parser.add_argument("--samples_per_target", type=int, default=100)
    parser.add_argument("--partial_samples_per_target", type=int, help="Partial flow: samples per target (default 8)")
    parser.add_argument("--framework_pdb", type=str, help="Framework PDB")
    parser.add_argument("--heavy_chain", type=str, help="Heavy chain id")
    parser.add_argument("--light_chain", type=str, help="Light chain id")
    parser.add_argument("--cdr_length", type=str, help="CDR length spec")

    # Sequence design knobs (paper-aligned defaults in config)
    parser.add_argument("--seq1_num_per_backbone", type=int, help="Seq1: sequences per backbone")
    parser.add_argument("--seq1_temp", type=float, help="Seq1: sampling temperature")
    parser.add_argument("--seq1_bias_large_residues", action="store_true", help="Seq1: bias bulky residues")
    parser.add_argument("--seq1_bias_num", type=int, help="Seq1: number of biased sequences")
    parser.add_argument("--seq1_bias_residues", type=str, help="Seq1: bias residues (comma-separated)")
    parser.add_argument("--seq1_bias_weight", type=float, help="Seq1: bias weight for residues")
    parser.add_argument("--seq2_num_per_backbone", type=int, help="Seq2: sequences per backbone")
    parser.add_argument("--seq2_temp", type=float, help="Seq2: sampling temperature")
    parser.add_argument("--seq2_use_soluble_ckpt", action="store_true", help="Seq2: use soluble checkpoint")

    # Antibody/VHH aliases
    parser.add_argument("--vhh_backbones", type=int, help="Alias for samples_per_target (antibody/vhh)")
    parser.add_argument("--vhh_cdr1_num", type=int, help="Alias for seq1_num_per_backbone (antibody/vhh)")
    parser.add_argument("--vhh_partial_num", type=int, help="Alias for partial_samples_per_target (antibody/vhh)")
    parser.add_argument("--vhh_cdr2_num", type=int, help="Alias for seq2_num_per_backbone (antibody/vhh)")

    # Filtering / ranking knobs
    parser.add_argument("--af3score1_iptm_min", type=float, help="AF3Score R1 ipTM threshold")
    parser.add_argument("--af3score1_ptm_min", type=float, help="AF3Score R1 pTM threshold")
    parser.add_argument("--af3score1_top_k", type=int, help="AF3Score R1: keep top K candidates (global)")
    parser.add_argument("--af3score2_iptm_min", type=float, help="AF3Score R2 ipTM threshold")
    parser.add_argument("--af3score2_ptm_min", type=float, help="AF3Score R2 pTM threshold")
    parser.add_argument("--af3score2_top_k", type=int, help="AF3Score R2: keep top K candidates (global)")
    parser.add_argument("--af3refold_iptm_min", type=float, help="AF3 refold ipTM threshold")
    parser.add_argument("--af3refold_ptm_min", type=float, help="AF3 refold pTM threshold")
    parser.add_argument("--af3refold_dockq_min", type=float, help="AF3 refold DockQ threshold")
    parser.add_argument("--af3refold_num_samples", type=int, help="AF3 refold: num_samples per seed")
    parser.add_argument("--af3refold_model_seeds", type=str, help="AF3 refold: model seeds (e.g. 0-19)")
    parser.add_argument(
        "--af3refold_no_templates",
        action="store_true",
        default=None,
        help="AF3 refold: disable templates",
    )
    parser.add_argument("--interface_energy_min", type=float, help="Rosetta interface energy cutoff (REU)")
    parser.add_argument("--interface_distance", type=float, help="Rosetta interface distance cutoff (A)")
    parser.add_argument("--relax_max_iter", type=int, help="Rosetta relax max iterations (legacy)")
    parser.add_argument("--relax_fixbb", action="store_true", default=None, help="Rosetta relax: fix backbone (legacy)")
    parser.add_argument(
        "--fixed_chains",
        type=str,
        help="Chains to fix backbone (underscore or comma separated, e.g. A_B or A,B)",
    )
    parser.add_argument("--dockq_min", type=float, help="DockQ minimum threshold")
    parser.add_argument("--partial_start_t", type=float, help="Partial flow start_t")
    parser.add_argument("--rank_top_k", type=int, help="Top K to keep in ranking")

    parser.add_argument("--af3_refold", action="store_true", help="Enable AF3 refold step")

    # Tool paths (optional)
    parser.add_argument("--ppiflow_ckpt", type=str, help="PPIFlow checkpoint path")
    parser.add_argument("--abmpnn_ckpt", type=str, help="AbMPNN checkpoint path")
    parser.add_argument("--mpnn_ckpt", type=str, help="ProteinMPNN checkpoint path")
    parser.add_argument("--mpnn_ckpt_soluble", type=str, help="ProteinMPNN soluble checkpoint path")
    parser.add_argument("--af3score_repo", type=str, help="AF3Score repo path")
    parser.add_argument("--rosetta_bin", type=str, help="Rosetta scripts binary path (rosetta_scripts)")
    parser.add_argument("--rosetta_db", type=str, help="Rosetta database path (optional)")
    parser.add_argument("--flowpacker_repo", type=str, help="FlowPacker repo path")
    parser.add_argument("--af3_weights", type=str, help="AF3 weights path")
    parser.add_argument("--mpnn_repo", type=str, help="ProteinMPNN repo path")
    parser.add_argument("--abmpnn_repo", type=str, help="AbMPNN repo path")
    parser.add_argument("--mpnn_run", type=str, help="Path to protein_mpnn_run.py")
    parser.add_argument("--abmpnn_run", type=str, help="Path to protein_mpnn_run.py for AbMPNN")
    parser.add_argument("--dockq_bin", type=str, help="Path to DockQ binary")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser("ppiflow")
    sub = parser.add_subparsers(dest="command", required=True)

    p_pipeline = sub.add_parser("pipeline", help="Configure then execute")
    _add_common_args(p_pipeline)
    _add_cli_input_args(p_pipeline)

    p_configure = sub.add_parser("configure", help="Write step configs and steps.yaml")
    _add_common_args(p_configure)
    _add_cli_input_args(p_configure)

    p_execute = sub.add_parser("execute", help="Execute from steps.yaml")
    p_execute.add_argument("--output", type=str, required=True)
    p_execute.add_argument("--steps", type=str, default="all")
    p_execute.add_argument("--reuse", action="store_true")
    p_execute.add_argument("--continue-on-error", action="store_true")
    p_execute.add_argument("--verbose", action="store_true")

    p_rank = sub.add_parser("rank", help="Run rank step only")
    p_rank.add_argument("--output", type=str, required=True)
    p_rank.add_argument("--reuse", action="store_true")
    p_rank.add_argument("--continue-on-error", action="store_true")
    p_rank.add_argument("--verbose", action="store_true")

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "pipeline":
        configure_pipeline(args)
        execute_pipeline(args)
    elif args.command == "configure":
        configure_pipeline(args)
    elif args.command == "execute":
        execute_pipeline(args)
    elif args.command == "rank":
        args.steps = "rank"
        execute_pipeline(args)
    else:
        raise SystemExit(f"Unknown command: {args.command}")


if __name__ == "__main__":
    main()
