import os
import glob
import argparse
import math
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--base_path', required=True)
parser.add_argument('--output_folder_base', required=True)
parser.add_argument('--num_of_jobs', type=int, required=True)
args = parser.parse_args()

pdbs = sorted(glob.glob(os.path.join(args.base_path, "*.pdb")))
N = len(pdbs)
n = args.num_of_jobs

batches = np.array_split(pdbs, n)

for idx, batch in enumerate(batches, 1):
    folder = os.path.join(args.output_folder_base, f"batch_{idx}")
    os.makedirs(folder, exist_ok=True)
    print(f"→ Batch {idx}: {len(batch)} files")
    for pdb in batch:
        dst = os.path.join(folder, os.path.basename(pdb))
        if not os.path.exists(dst):
            os.link(pdb, dst)
