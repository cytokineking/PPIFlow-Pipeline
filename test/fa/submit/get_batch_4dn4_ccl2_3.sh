#!/bin/bash 
#SBATCH -p gpu41,gpu43 
#SBATCH --gres=gpu:1 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1
#SBATCH -N 1 
#SBATCH -o /lustre/grp/cmclab/qinxy/antibody_ppiflow/demo/test/fa/logs/%j-4dn4_ccl2_3.out
#SBATCH -J 3_proteinmpnn

echo "Job started at: $(date)"
start_time=$(date +%s)
export WANDB_MODE=disabled
source /lustre/grp/cmclab/qinxy/new_miniconda3/etc/profile.d/conda.sh 
conda activate mlfold 
folder_with_pdbs=/lustre/grp/cmclab/qinxy/antibody_ppiflow/demo/test/pf/4dn4_ccl2_3 
fixed_positions="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 115 116 117 118 119 120 121 122 123 124 125" 
output_dir=/lustre/grp/cmclab/qinxy/antibody_ppiflow/demo/test/fa/4dn4_ccl2_3 

path_for_parsed_chains=/lustre/grp/cmclab/qinxy/antibody_ppiflow/demo/test/fa/4dn4_ccl2_3/parsed_pdbs.jsonl 
path_for_assigned_chains=/lustre/grp/cmclab/qinxy/antibody_ppiflow/demo/test/fa/4dn4_ccl2_3/assigned_pdbs.jsonl 

path_for_fixed_positions=/lustre/grp/cmclab/qinxy/antibody_ppiflow/demo/test/fa/4dn4_ccl2_3/fixed_pdbs.jsonl 

python ~/software/ProteinMPNN-main/helper_scripts/parse_multiple_chains.py --input_path=/lustre/grp/cmclab/qinxy/antibody_ppiflow/demo/test/pf/4dn4_ccl2_3 --output_path=$path_for_parsed_chains 
python ~/software/ProteinMPNN-main/helper_scripts/assign_fixed_chains.py --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list A 

python ~/software/ProteinMPNN-main/helper_scripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list A --position_list '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 115 116 117 118 119 120 121 122 123 124 125' 
python ~/software/ProteinMPNN-main/protein_mpnn_run.py \
--path_to_model_weights /lustre/grp/cmclab/share/wanghz/Programs/ProteinMPNN-qin/model_weights \
--model_name 'abmpnn' \
--jsonl_path $path_for_parsed_chains \
--chain_id_jsonl $path_for_assigned_chains \
--fixed_positions_jsonl $path_for_fixed_positions \
--out_folder $output_dir \
--num_seq_per_target 4 \
--sampling_temp 0.1 \
--batch_size 4 \
--use_soluble_model \
--omit_AAs C 

end_time=$(date +%s)
echo "Job ended at: $(date)"
elapsed=$((end_time - start_time))
echo "Elapsed time: $elapsed seconds"