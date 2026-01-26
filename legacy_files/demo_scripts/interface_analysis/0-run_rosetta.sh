#!/bin/bash
set -euo pipefail

# Usage check
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir>"
    exit 1
fi

input_dir="$1"
output_dir="$2"

# Resolve absolute paths
input_dir="$(readlink -f "$input_dir")"
output_dir="$(readlink -f "$output_dir")"

# Define the script working directory (where codes/ reside)
# Adjust this if the python/sh scripts are located elsewhere
WD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$output_dir"

# --- Enumerate PDB files ----------------------------------------------
mapfile -t pdb_files < <(find "$input_dir" -maxdepth 1 -name '*.pdb' | sort)
echo "Found ${#pdb_files[@]} PDB files."

# --- Main loop --------------------------------------------------------
for pdb_path in "${pdb_files[@]}"; do
    folder_name="$(basename "$pdb_path" .pdb)"
    job_dir="$output_dir/$folder_name"

    # Skip if the output directory already exists
    if [ -d "$job_dir" ]; then
        echo "Skipping $folder_name: Output directory already exists (assumed complete)."
        continue
    fi

    echo -e "\nCreating task: $folder_name"
    mkdir -p "$job_dir"/{out,plots}
    cp "$pdb_path" "$job_dir/${folder_name}.pdb"

    # Enter job directory
    pushd "$job_dir" >/dev/null

    # Absolute path variables
    pdb_abs="$job_dir/${folder_name}.pdb"
    out_path="$job_dir/out/${folder_name}.out"
    log_path="$job_dir/out/${folder_name}.log"

    echo "Working Directory: $(pwd)"
    echo "PDB Path:          $pdb_abs"
    echo "Log Path:          $log_path"

    # Run XML update script
    /lustre/grp/cmclab/qinxy/new_miniconda3/envs/py39/bin/python "$WD/codes/update_xml.py" --pdb_path "$pdb_abs"

    # Submit SLURM job
    sbatch "$WD/codes/per-residue-energy_submit.sh" \
        "update.xml" "$pdb_abs" "$out_path"

    popd >/dev/null
done

echo -e "\nAll tasks have been successfully prepared in $output_dir."
