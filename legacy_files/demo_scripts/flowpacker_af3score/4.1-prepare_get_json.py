from Bio import PDB
from Bio.PDB import Structure, Model
from Bio.PDB import PDBParser, MMCIFIO
import os
import pandas as pd
from tqdm import tqdm
import multiprocessing as mp
from pathlib import Path
import math
import shutil
import argparse
import multiprocessing as mp
from concurrent.futures import ThreadPoolExecutor, as_completed
import json


def format_msa_sequence(sequence):
    return f">query\n{sequence}\n"

def get_chain_sequences(row):
    chain_sequences = []
    chain_columns = [col for col in row.index if col.startswith('chain_') and col.endswith('_seq')]
    for col in chain_columns:
        if pd.notna(row[col]) and row[col] != '':
            chain_id = col.split('_')[1]
            chain_sequences.append((chain_id, row[col]))
    return chain_sequences



protein_letters_3to1 = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
    'MSE': 'M'
}

def get_sequence(chain):
    sequence = ""
    for residue in chain:
        if residue.id[0] == ' ':
            resname = residue.get_resname().upper()
            sequence += protein_letters_3to1.get(resname, 'X')
    return sequence

def process_single_pdb(args):
    input_pdb, output_dir_cif = args
    try:
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("structure", input_pdb)
        base_name = os.path.splitext(os.path.basename(input_pdb))[0]

        chain_sequences = {}
        merged_sequence = ""

        for chain in structure[0]:
            chain_id = chain.id
            sequence = get_sequence(chain)
            chain_sequences[chain_id] = sequence
            merged_sequence += sequence

            new_structure = Structure.Structure("new_structure")
            new_model = Model.Model(0)
            new_structure.add(new_model)
            new_model.add(chain.copy())

            cif_io = MMCIFIO()
            cif_io.set_structure(new_structure)
            cif_output = os.path.join(output_dir_cif, f"{base_name}_chain_{chain_id}.cif")
            cif_io.save(cif_output)

        return base_name, chain_sequences, len(merged_sequence)
    
    except Exception as e:
        print(f"Error processing {input_pdb}: {str(e)}")
        return None, None, None

def generate_json_files(tasks):
    row, cif_dir, output_dir = tasks
    complex_name = row['complex']
    chain_sequences = get_chain_sequences(row)
    
    if not chain_sequences:
        print(f"⚠️ 警告: {complex_name} 没有有效的chain序列")
        return None

    sequences = []
    for chain_id, sequence in chain_sequences:
        cif_filename = f"{complex_name}_chain_{chain_id}.cif"
        cif_path = os.path.join(cif_dir, cif_filename)
        if not os.path.exists(cif_path):
            print(f"⚠️ 警告: {cif_filename} 不存在")
            continue
        sequences.append({
            "protein": {
                "id": chain_id,
                "sequence": sequence,
                "modifications": [],
                "unpairedMsa": format_msa_sequence(sequence),
                "pairedMsa": format_msa_sequence(sequence),
                "templates": [{
                    "mmcifPath": cif_path,
                    "queryIndices": list(range(len(sequence))),
                    "templateIndices": list(range(len(sequence)))
                }]
            }
        })
    
    if not sequences:
        print(f"⚠️ 警告: {complex_name} 没有有效的序列数据")
        return None
    
    json_data = {
        "dialect": "alphafold3",
        "version": 1,
        "name": complex_name,
        "sequences": sequences,
        "modelSeeds": [10],
        "bondedAtomPairs": None,
        "userCCD": None
    }
    
    output_filename = f"{complex_name}.json"
    output_path = os.path.join(output_dir, output_filename)
    with open(output_path, 'w') as f:
        json.dump(json_data, f, indent=2)
    
    chain_ids = [chain[0] for chain in chain_sequences]
    # print(f"✅ 生成JSON文件: {output_filename} (链: {', '.join(chain_ids)})")
    return output_filename

    # 统计成功生成的数量
    success_count = sum(1 for r in results if r is not None)
    print(f"\n✅ 生成完毕，共生成 {success_count} 个JSON文件")

    
def get_seq_main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, default='/lustre/grp/cmclab/share/wuxy/af3score/scripts/merge_process/test_data/input_pdbs')
    parser.add_argument('--output_dir_cif', type=str, default='/lustre/grp/cmclab/share/wuxy/af3score/scripts/merge_process/test_data/single_chain_cif')
    parser.add_argument('--save_csv', type=str, default='/lustre/grp/cmclab/share/wuxy/af3score/scripts/merge_process/test_data/seq.csv')
    parser.add_argument('--output_dir_json', type=str, default='/lustre/grp/cmclab/share/wuxy/af3score/scripts/merge_process/test_data/json')
    parser.add_argument('--batch_dir', type=str, default='/lustre/grp/cmclab/share/wuxy/process/flowpacker_af3score_v3/data/test_1/af3score_base_outputs/batch')
    parser.add_argument('--num_workers', type=int, default=36, help='使用的进程数（默认36）')
    parser.add_argument('--num_jobs', type=int, default=None)
    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir_cif = args.output_dir_cif
    output_dir_json = args.output_dir_json
    save_csv = args.save_csv
    batch_dir = args.batch_dir
    num_workers = args.num_workers
    num_jobs = args.num_jobs

    print(f'input_dir: {input_dir}')
    print(f'output_dir_cif: {output_dir_cif}')
    print(f'output_dir_json: {output_dir_json}')

    os.makedirs(output_dir_cif, exist_ok=True)

    pdb_files = list(Path(input_dir).glob("*.pdb"))
    args = [(str(f), output_dir_cif) for f in pdb_files]

    sequences_dict = {}
    with mp.Pool(processes=mp.cpu_count()) as pool:
        results = list(tqdm(pool.imap(process_single_pdb, args), total=len(pdb_files), desc="Processing PDB files"))

    for base_name, chain_sequences, length in results:
        if base_name is not None:
            sequences_dict[base_name] = {
                "sequences": chain_sequences,
                "length": length
            }

    all_chain_ids = set()
    for entry in sequences_dict.values():
        all_chain_ids.update(entry["sequences"].keys())


    def chain_sort_key(chain_id):
        return str(chain_id)

    all_chain_ids = sorted(list(all_chain_ids), key=chain_sort_key)

    rows = []
    for complex_name, entry in sequences_dict.items():
        chain_data = entry["sequences"]
        length = entry["length"]
        row = {'complex': complex_name, 'total_length': length}
        for chain_id in all_chain_ids:
            row[f'chain_{chain_id}_seq'] = chain_data.get(chain_id, '')
        rows.append(row)

    df = pd.DataFrame(rows)
    cols = ['complex', 'total_length'] + [col for col in df.columns if col not in ['complex', 'total_length']]
    df = df[cols]

    df.to_csv(save_csv, index=False)
    print(f"\n✅ 序列信息已保存到 {save_csv}")

    rows = [row for _, row in df.iterrows()]
    tasks1 = [(r, output_dir_cif, output_dir_json) for _, r in df.iterrows()]
    with mp.Pool(processes=num_workers) as pool:
        results = list(tqdm(pool.imap(generate_json_files, tasks1), total=len(tasks1), desc="Processing Rows"))

    success_count = sum(1 for r in results if r)
    print(f"\n✅ 并行生成JSON文件完成，共生成：{success_count} 个文件")


    batch_json = f'{batch_dir}/json'
    batch_pdb = f'{batch_dir}/pdb'

    df = df.sort_values("total_length").reset_index(drop=True)
    bz = math.ceil(len(df) / num_jobs)
    print(f"📦 总共 {len(df)} 个样本，分 {num_jobs} 个 batch，每个 batch 最多 {bz} 个")

    for i in range(num_jobs):
        sub = df.iloc[i * bz:(i + 1) * bz]
        if sub.empty:
            print(f"🔹 Batch {i} 为空，跳过")
            continue

        mx = sub['total_length'].max()
        name = f"batch_{i}_{mx}"
        bd_json = os.path.join(batch_json, name)
        bd_pdb = os.path.join(batch_pdb, name)
        os.makedirs(bd_json, exist_ok=True)
        os.makedirs(bd_pdb, exist_ok=True)

        count = 0 
        for _, r in sub.iterrows():
            cid = r['complex']
            for ext, od in [('.pdb', input_dir), ('.json', output_dir_json)]:
                src = os.path.join(od, f"{cid}{ext}")
                if os.path.exists(src) and src.endswith('.json'):
                    shutil.copy(src, bd_json)
                else:
                    shutil.copy(src, bd_pdb)
            count += 1

        print(f"✅ {name}: 包含 {count} 个 complex")

    print(f"\n📊 总共处理 {len(df)} 个 complex，生成了 {len(os.listdir(batch_json))} 个 batch")

if __name__ == "__main__":
    get_seq_main()

