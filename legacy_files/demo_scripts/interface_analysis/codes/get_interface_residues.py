import numpy as np
import os
import warnings
from Bio import PDB
import argparse

# Suppress specific warnings
warnings.filterwarnings("ignore")

def get_interface_res(pdb, dist_cutoff, chain1, chain2):
    parser = PDB.PDBParser()

    pdb_name = os.path.basename(pdb).replace('.pdb', '')
    structure = parser.get_structure(pdb_name, pdb)
    model = structure[0]
    chain_A = model[chain1]
    chain_B = model[chain2]
    
    def get_coords(chain):
        coords = []
        res_ids = []
        for res in chain:
            # Try to get 'CA', if not present, try 'CB'
            atom = res['CA'] if res.has_id('CA') else (res['CB'] if res.has_id('CB') else None)
            if atom:
                coords.append(atom.get_coord())
                res_ids.append(res.get_id()[1])
            else:
                print(f"Warning: Residue {res.get_id()} in chain {chain.get_id()} does not have 'CA' or 'CB' atom.")
        return np.array(coords).T, res_ids

    chain_A_CA, chain_A_res_ids = get_coords(chain_A)
    chain_B_CA, chain_B_res_ids = get_coords(chain_B)
    
    dist = np.sqrt(np.sum((chain_A_CA[:, :, None] - chain_B_CA[:, None, :])**2, axis=0))
    
    interface_residues = np.where(dist < dist_cutoff)

    # return interface_residues
    return [chain_A_res_ids[i] for i in interface_residues[0]], [chain_B_res_ids[i] for i in interface_residues[1]]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process find residues between interface.")
    parser.add_argument('-i', '--input', required=True, help='Path to the main folder containing subdirectories with PDB files and PLUMED input files.')
    parser.add_argument('-c1', '--chain1', required=True, help='ID of the first chain.')
    parser.add_argument('-c2', '--chain2', required=True, help='ID of the second chain.')
    parser.add_argument('-d', '--cutoff', type=float, default=10.0, help='Cutoff distance for determining interface atoms (default: 10.0).')

    args = parser.parse_args()
    interface_residues = get_interface_res(args.input, args.cutoff, args.chain1, args.chain2)
    unique_residues_A = set(interface_residues[0])
    unique_residues_B = set(interface_residues[1])
    
    # Print residues in the desired format
    print(f"Interface residues in chain {args.chain1}: {','.join(map(str, sorted(unique_residues_A)))}")
    print(f"Interface residues in chain {args.chain2}: {','.join(map(str, sorted(unique_residues_B)))}")
