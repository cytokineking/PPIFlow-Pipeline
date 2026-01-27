#!/usr/bin/env python
"""
Wrapper around the OG get_interface_energy.py that patches the pandas view/copy
issue so binder_key/target_key are available on interface_score_df.
This keeps OG logic intact while avoiding empty outputs under newer pandas.
"""
from __future__ import annotations

import argparse
import os
from functools import partial
from multiprocessing import Pool, cpu_count
from pathlib import Path
import importlib.util


def _load_og_module(path: Path):
    spec = importlib.util.spec_from_file_location("ppiflow_og_get_interface_energy", str(path))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Failed to load OG script at {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _patched_get_interface_energy(og):
    def _impl(interchain_score_df, interface_pair, binder_id, target_id, plot_path):
        if interchain_score_df is None or len(interchain_score_df) == 0:
            return {}

        binder_ids = og._parse_chain_ids(binder_id)
        target_ids = og._parse_chain_ids(target_id)

        df = interchain_score_df
        if binder_ids:
            df = df.loc[df['binder_id'].isin(binder_ids)]
        if target_ids:
            df = df.loc[df['target_id'].isin(target_ids)]

        df = df.copy()
        if len(binder_ids) == 1:
            df['binder_key'] = df['binder_res'].astype(int)
        else:
            df['binder_key'] = df['binder_id'].astype(str) + df['binder_res'].astype(str)
        if len(target_ids) == 1:
            df['target_key'] = df['target_res'].astype(int)
        else:
            df['target_key'] = df['target_id'].astype(str) + df['target_res'].astype(str)

        df['in_interface'] = df.apply(
            lambda row: (row['binder_id'], int(row['binder_res']), row['target_id'], int(row['target_res'])) in interface_pair,
            axis=1,
        )
        interface_score_df = df.loc[df['in_interface'] == True]

        summed_df = interface_score_df.groupby('binder_key')['total'].sum().reset_index()
        summed_dict = summed_df.set_index('binder_key')['total'].to_dict()
        print(interface_score_df.head(2))
        print(interface_score_df.columns)
        og.plot_score(df, plot_path)
        return summed_dict

    return _impl


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_csv', type=str, help='csv file contains pdb pdbpath. When the input_pdbdir is not provided, input_csv must be provided!')
    parser.add_argument('--input_pdbdir', type=str, help='The pdb folder contains pdb files. When the input_csv is not provided, input_pdbdir must be provided!')
    parser.add_argument('--rosetta_dir', type=str, help='rosetta results directory of all pdb results are saved', required=True)
    parser.add_argument('--binder_id', type=str, help='binder chain id', required=True)
    parser.add_argument('--target_id', type=str, help='target chain id', required=True)
    parser.add_argument('--output_dir', type=str, help='output csv file of all pdb results', required=True)
    parser.add_argument("--interface_dist", type=float, default=12.0, help="interface distance between target and binder")
    parser.add_argument("--og_script", type=str, default="", help="Path to OG get_interface_energy.py")
    args = parser.parse_args()

    if args.og_script:
        og_path = Path(args.og_script).resolve()
    else:
        og_path = Path(__file__).resolve().parents[1] / "PPIFlow_OG_repo" / "demo_scripts" / "interface_analysis" / "get_interface_energy.py"
    if not og_path.exists():
        fallback = Path(__file__).resolve().parents[1] / "legacy" / "demo_scripts" / "interface_analysis" / "get_interface_energy.py"
        if fallback.exists():
            og_path = fallback
        else:
            raise SystemExit(f"OG script not found at {og_path}")

    og = _load_og_module(og_path)
    og.get_interface_energy = _patched_get_interface_energy(og)

    output_dir = args.output_dir
    interface_dist = args.interface_dist
    os.makedirs(os.path.dirname(output_dir), exist_ok=True)

    df = og.get_input_df(args)

    print("start serial")
    func = partial(og.main, output_dir=output_dir, distance_threshold=interface_dist)
    results = [func(row) for _, row in og.tqdm(df.iterrows(), total=len(df))]

    df['binder_energy'] = results
    df.drop(columns=['pdbpath', 'pdbname'])
    df.to_csv(os.path.join(output_dir, "residue_energy.csv"), index=False)
    print(output_dir)

    title = f"interface_binder_residues_energy_sum"
    savepath = os.path.join(output_dir, 'residue_energy_interface_binder_residues_energy_sum.png')
    og.plot_binder_score(os.path.join(output_dir, "residue_energy.csv"), title=title, savepath=savepath)


if __name__ == '__main__':
    main()
