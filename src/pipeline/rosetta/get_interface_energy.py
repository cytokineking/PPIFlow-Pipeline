import os
from glob import glob
from functools import partial
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from Bio import PDB
import ast
import argparse
import warnings
warnings.filterwarnings("ignore")


def get_rosetta_result(logfile):
    outlist = []
    with open(logfile, 'r') as infile:
        for line in infile:
            if line.rstrip('\n').startswith("ResResE"):
                outlist.append(line.split())
    if len(outlist) == 0:
        print(logfile)
        print("No logfile Error")
        return None
    outdf = pd.DataFrame(outlist[2:], columns=outlist[0])
    return outdf


def get_interchain_score(rosetta_path):
    raw_score_df = get_rosetta_result(rosetta_path)
    if raw_score_df is None:
        return None

    chain_map = {str(i): chr(i) for i in list(range(65, 91)) + list(range(97, 123))}
    raw_score_df['binder_id'] = raw_score_df['Res1'].str.split("_", expand=True)[1].str[:2].map(chain_map)
    raw_score_df['target_id'] = raw_score_df['Res2'].str.split("_", expand=True)[1].str[:2].map(chain_map)
    raw_score_df['binder_res'] = raw_score_df['Res1'].str.split("_", expand=True)[1].str[2:]
    raw_score_df['target_res'] = raw_score_df['Res2'].str.split("_", expand=True)[1].str[2:]

    interchain_score_df = raw_score_df.loc[raw_score_df['binder_id'] != raw_score_df['target_id'],
            ['binder_id', 'target_id', 'binder_res', 'target_res', 'total']]
    interchain_score_df['total'] = interchain_score_df['total'].astype(float)
    interchain_score_df = interchain_score_df.loc[interchain_score_df['total'] < 0]
    interchain_score_df.reset_index(drop=True, inplace=True)
    return interchain_score_df


def get_cb_or_ca(residue):
    if "CB" in residue:
        return residue["CB"].coord
    if "CA" in residue:
        return residue["CA"].coord
    return None


def _parse_chain_ids(value):
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return [str(v).strip() for v in value if str(v).strip()]
    if not isinstance(value, str):
        value = str(value)
    value = value.replace("_", ",")
    return [v.strip() for v in value.split(",") if v.strip()]


def get_residue_pairs_within_distance(pdb_file, binder_id, target_id, distance_threshold=10.0):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    model = structure[0]
    binder_ids = _parse_chain_ids(binder_id)
    target_ids = _parse_chain_ids(target_id)

    selected_pairs = set()
    for b_chain in binder_ids:
        if not model.has_id(b_chain):
            warnings.warn(f"Binder chain {b_chain} not found in {pdb_file}")
            continue
        binder_chain = model[b_chain]
        for t_chain in target_ids:
            if not model.has_id(t_chain):
                warnings.warn(f"Target chain {t_chain} not found in {pdb_file}")
                continue
            target_chain = model[t_chain]
            for res1 in binder_chain:
                coord1 = get_cb_or_ca(res1)
                for res2 in target_chain:
                    coord2 = get_cb_or_ca(res2)
                    if coord1 is not None and coord2 is not None:
                        distance = np.linalg.norm(coord1 - coord2)
                        if distance <= distance_threshold:
                            selected_pairs.add((b_chain, res1.id[1], t_chain, res2.id[1]))
    return selected_pairs


def plot_score(df, plot_path):
    if "target_key" in df.columns and "binder_key" in df.columns:
        heatmap_data = df.pivot(index="target_key", columns="binder_key", values="total")
    else:
        heatmap_data = df.pivot(index="target_res", columns="binder_res", values="total")

    plt.figure(figsize=(6, 5))
    sns.heatmap(heatmap_data, cmap="coolwarm", fmt=".1f", vmin=-5, vmax=1)

    plt.title("Score Heatmap")
    plt.xlabel("binder")
    plt.ylabel("target")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=100)


def get_interface_energy(interchain_score_df, interface_pair, binder_id, target_id, plot_path):
    if interchain_score_df is None or len(interchain_score_df) == 0:
        return {}

    binder_ids = _parse_chain_ids(binder_id)
    target_ids = _parse_chain_ids(target_id)

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
    plot_score(df, plot_path)
    return summed_dict


def get_input_df(args):
    input_csv = args.input_csv
    input_pdbdir = args.input_pdbdir
    rosetta_dir = f"{args.rosetta_dir}"

    if (input_csv is None) and (input_pdbdir is None):
        raise ValueError("input_csv and input_pdbdir cannot be None at the same time")
    if (input_csv is not None) and (input_pdbdir is not None):
        raise ValueError("Only one input_csv and input_pdbdir needs to be provided.")
    if input_csv:
        df = pd.read_csv(input_csv)
    else:
        pdbfiles = glob(f"{input_pdbdir}/*.pdb")
        print(f"Found pdb: {len(pdbfiles)}")
        df = pd.DataFrame({'pdbpath': pdbfiles})
    df["pdbname"] = df["pdbpath"].apply(lambda x: os.path.splitext(os.path.basename(x))[0])
    df['rosetta_path'] = df['pdbpath'].apply(lambda
                                              x: f"{rosetta_dir}/{os.path.basename(x).split('.pdb')[0]}/out/{os.path.basename(x).split('.pdb')[0]}.out")

    df["target_id"] = args.target_id
    df['binder_id'] = args.binder_id

    return df


def main(row, output_dir="", distance_threshold=10):
    logfile = row['rosetta_path']
    pdb_file = row['pdbpath']
    binder_id = row['binder_id']
    target_id = row['target_id']
    try:
        if not os.path.exists(logfile):
            print(row)
            print('logfile does not exist')
        interchain_score_path = get_interchain_score(logfile)
        interface_pair = get_residue_pairs_within_distance(pdb_file, binder_id, target_id, distance_threshold=distance_threshold)
        print(os.path.join(output_dir, os.path.basename(pdb_file).split('.pdb')[0] + ".png"))
        summed_dict = get_interface_energy(
            interchain_score_path,
            interface_pair,
            binder_id,
            target_id,
            os.path.join(output_dir, os.path.basename(pdb_file).split('.pdb')[0] + ".png"),
        )
    except Exception as e:
        print(e)
        return {}

    return summed_dict


def plot_binder_score(df_path, title, fontsize=15, savepath=None):
    df = pd.read_csv(df_path)
    df['binder_energy'] = df['binder_energy'].apply(lambda x: ast.literal_eval(x))
    binder_energy = [value for dictionary in df['binder_energy'] for value in dictionary.values()]
    print("Extracted values:", len(binder_energy))
    plt.figure(figsize=(10, 6))
    sns.kdeplot(binder_energy, fill=True)

    plt.title(title, fontsize=15)
    plt.tick_params(axis='x', labelsize=fontsize)
    plt.tick_params(axis='y', labelsize=fontsize)
    plt.xlabel('energy', fontsize=fontsize)
    plt.ylabel('density', fontsize=fontsize)

    plt.xlim((-10, 0))
    plt.tight_layout()
    if savepath is not None:
        plt.savefig(savepath)

    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_csv', type=str, help='csv file contains pdb pdbpath. When the input_pdbdir is not provided, input_csv must be provided!')
    parser.add_argument('--input_pdbdir', type=str, help='The pdb folder contains pdb files. When the input_csv is not provided, input_pdbdir must be provided!')
    parser.add_argument('--rosetta_dir', type=str, help='rosetta results directory of all pdb results are saved', required=True)
    parser.add_argument('--binder_id', type=str, help='binder chain id', required=True)
    parser.add_argument('--target_id', type=str, help='target chain id', required=True)
    parser.add_argument('--output_dir', type=str, help='output csv file of all pdb results', required=True)
    parser.add_argument("--interface_dist", type=float, default=12.0, help="interface distance between target and binder")

    args = parser.parse_args()

    output_dir = args.output_dir
    interface_dist = args.interface_dist
    os.makedirs(os.path.dirname(output_dir), exist_ok=True)

    df = get_input_df(args)

    print("start Pool")
    func = partial(main, output_dir=output_dir, distance_threshold=interface_dist)
    with Pool(processes=cpu_count()) as pool:
        results = list(tqdm(pool.imap(func, [row for _, row in df.iterrows()]), total=len(df)))

    df['binder_energy'] = results

    df.drop(columns=['pdbpath', 'pdbname'])
    df.to_csv(os.path.join(output_dir, "residue_energy.csv"), index=False)
    print(output_dir)

    title = "interface_binder_residues_energy_sum"
    savepath = os.path.join(output_dir, 'residue_energy_interface_binder_residues_energy_sum.png')
    plot_binder_score(os.path.join(output_dir, "residue_energy.csv"), title=title, savepath=savepath)
