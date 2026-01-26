import os,sys
import csv
from multiprocessing import Pool

three_to_one = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def extract_sequence_from_pdb(pdb_file):
    sequences = {}
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") and line[13:15].strip() == "CA":
                chain_id = line[21].strip()
                res_name = line[17:20].strip()
                if res_name in three_to_one:
                    if chain_id not in sequences:
                        sequences[chain_id] = []
                    sequences[chain_id].append(three_to_one[res_name])

    sequence_str = ':'.join([''.join(seq) for seq in sequences.values()])
    return os.path.basename(pdb_file), sequence_str

def write_fa(result, output_path):
    with open(output_path, 'w') as fa_file:
        for pdb_name, sequence in result:
            seq_name = pdb_name.rsplit('.', 1)[0]
            if sequence:
                fa_file.write(f">{seq_name}\n{sequence}\n")

def write_csv(result, output_path):
    with open(output_path, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['pdb_name', 'sequence'])
        for pdb_name, sequence in result:
            if sequence:
                writer.writerow([pdb_name, sequence])

def process_pdb_files(folder_path, output_path, output_format='fa', num_processes=4):
    pdb_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith(".pdb")]

    with Pool(processes=num_processes) as pool:
        result = pool.map(extract_sequence_from_pdb, pdb_files)

    if output_format == 'fa':
        write_fa(result, output_path)
    elif output_format == 'csv':
        write_csv(result, output_path)
    else:
        raise ValueError("Unsupported output format. Use 'fa' or 'csv'.")

if __name__ == "__main__":
    folder_path = sys.argv[1]
    output_format = sys.argv[2]  # 'fa' or 'csv'
    output_path = sys.argv[3]
    # output_path = f"{folder_path}_all.{output_format}"
    num_processes = 64

    process_pdb_files(folder_path, output_path, output_format, num_processes)
    print(f'seq saved at {output_path}')