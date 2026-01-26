#!/bin/bash
#SBATCH -p cpu1,cpu2,cpu3,fat,hygon
#SBATCH -n 1
#SBATCH -o plot-%j.out
#SBATCH -J plot

python /lustre/grp/cmclab/qinxy/script/interface_analysis/codes/plot_heatmap.py \
--log_path $1 --pdb_name $2 --plot_path $3 --chain1 $4 --chain2 $5 --residues1 $6 --residues2 $7

# log_path = args.log_path
# pdb_name = args.pdb_name
# plot_path = args.plot_path
# chain1 = args.chain1
# chain2 = args.chain2
# residues1 = args.residues1
# residues2 = args.residues2