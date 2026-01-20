
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

    
def get_seq_main():

    
    parser = argparse.ArgumentParser()
    # parser.add_argument('config', type=str)
    parser.add_argument('--input_dir', type=str, default='/lustre/grp/cmclab/share/wuxy/af3score/scripts/merge_process/test_data/input_pdbs')
    parser.add_argument('--output_dir_cif', type=str, default='/lustre/grp/cmclab/share/wuxy/af3score/scripts/merge_process/test_data/single_chain_cif')
    parser.add_argument('--save_csv', type=str, default='/lustre/grp/cmclab/share/wuxy/af3score/scripts/merge_process/test_data/seq.csv')
    parser.add_argument('--output_dir_json', type=str, default='/lustre/grp/cmclab/share/wuxy/af3score/scripts/merge_process/test_data/json')
    parser.add_argument('--batch_dir', type=str, default='/lustre/grp/cmclab/share/wuxy/process/flowpacker_af3score_v3/data/test_1/af3score_base_outputs/batch')
    parser.add_argument('--num_workers', type=int, default=36)
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

    batch_json = f'{batch_dir}/json'
    batch_pdb = f'{batch_dir}/pdb'
    df=pd.read_csv('round_2/20250912_122246/af3score_base_outputs/single_seq.csv')
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
            for ext, od, bd in [('.pdb', input_dir, bd_pdb), ('.json', output_dir_json, bd_json)]:
                src = os.path.join(od, f"{cid}{ext}")
                dst = os.path.join(bd, f"{cid}{ext}") 
                try:
                    if not os.path.exists(dst):
                        os.symlink(src, dst)
                except FileExistsError:
                    pass
                except Exception as e:
                    print(f"⚠️ link 失败: {src} -> {dst}, 错误: {e}")


        print(f"✅ {name}: 包含 {count} 个 complex")

    print(f"\n📊 总共处理 {len(df)} 个 complex，生成了 {len(os.listdir(batch_json))} 个 batch")

if __name__ == "__main__":
    get_seq_main()

