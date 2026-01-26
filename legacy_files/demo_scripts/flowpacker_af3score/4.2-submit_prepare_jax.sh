#!/bin/bash
#SBATCH -p gpu41
#SBATCH --gres=gpu:1
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -J jax



/lustre/grp/cmclab/share/anaconda3/envs/jaxaf2/bin/python demo_scripts/flowpacker_af3score/4.2_prepare_pdb2jax.py \
      --pdb_folder $1 \
      --output_folder $2 
