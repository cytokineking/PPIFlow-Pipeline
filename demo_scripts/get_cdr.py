import os,sys
import pandas as pd
from Bio import PDB
from Bio.SeqUtils import seq1
from abnumber import Chain, ChainParseError
import concurrent.futures
from tqdm import tqdm
import warnings

warnings.filterwarnings("ignore")

def get_seq_and_residues(chain_obj):

    residues = []
    seq_chars = []
    for res in chain_obj:
        if PDB.is_aa(res, standard=True):
            try:
                aa = seq1(res.get_resname())
                seq_chars.append(aa)
                residues.append(res)
            except Exception:
                continue
    return "".join(seq_chars), residues

def format_res_id(res):

    _, num, insert = res.get_id()
    return f"{num}{insert.strip()}"

def process_single_file(args):

    pdb_path, scheme, cdr_definition = args
    pdb_id = os.path.basename(pdb_path).replace('.pdb', '').replace('.cif', '')
    
    parser = PDB.PDBParser(QUIET=True)
    
    try:
        structure = parser.get_structure(pdb_id, pdb_path)
    except Exception:
        return None 

    result_entry = None

    for model in structure:

        if 'A' in model:
            bio_chain = model['A']
        else:
            break
            
        full_seq, bio_residues = get_seq_and_residues(bio_chain)
        
        if len(full_seq) < 50: 
            break

        try:
            ab_chain = Chain(full_seq, scheme=scheme, cdr_definition=cdr_definition)
        except (ChainParseError, Exception):
            break

        v_seq_str = ab_chain.seq
        start_index = full_seq.find(v_seq_str)
        if start_index == -1:
            break
            
        v_seq_positions = list(ab_chain.positions.keys())
        
        fw_indices = []
        cdr_indices = []
        

        regions_order = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4']
        
        for region_name in regions_order:

            region_dict = getattr(ab_chain, f"{region_name.lower()}_dict", None)
            
            if not region_dict:
                continue
            
            is_cdr = 'CDR' in region_name.upper()
            
            for pos_obj in region_dict.keys():
                try:
                    idx_in_v = v_seq_positions.index(pos_obj)
                    
                    abs_idx = start_index + idx_in_v
                    
                    res = bio_residues[abs_idx]
                    
                    res_str = format_res_id(res)
                    
                    if is_cdr:
                        cdr_indices.append(res_str)
                    else:
                        fw_indices.append(res_str)
                        
                except (ValueError, IndexError):
                    continue
        
        if fw_indices or cdr_indices:
            result_entry = {
                'pdb_name': pdb_id,
                'fw_index': " ".join(fw_indices),
                'r2_cdr_pos': ",".join(f"A{i}" for i in cdr_indices)
            }
        
        break
        
    return result_entry

def run_parallel(pdb_files, scheme='imgt', cdr_definition='imgt', max_workers=None):
    """
    并行调度主函数
    """
    tasks = [(f, scheme, cdr_definition) for f in pdb_files]
    final_data = []
    
    if max_workers is None:
        max_workers = max(1, os.cpu_count() - 2)
    
    print(f"Starting processing {len(pdb_files)} files with {max_workers} cores...")
    print("Extracting Framework and CDR indices for Chain A...")
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        results_generator = list(tqdm(executor.map(process_single_file, tasks), total=len(tasks), unit="pdb"))
        
        for res in results_generator:
            if res: 
                final_data.append(res)
            
    return pd.DataFrame(final_data)


if __name__ == '__main__':
    
    import glob

    pdb_folder = sys.argv[1]  
    pdb_files = glob.glob(os.path.join(pdb_folder, "*.pdb"))
    
    if not pdb_files:
        print(f"No PDB files were found in {pdb_folder}.")
    else:

        df_final = run_parallel(pdb_files, scheme='imgt', cdr_definition='imgt')
        

        output_csv = sys.argv[2]  
        
        if not df_final.empty:
            df_final = df_final[['pdb_name', 'fw_index', 'r2_cdr_pos']]
            df_final.to_csv(output_csv, index=False)
            print(f"\nSaved at {output_csv}")
        else:
            print("ERROR")