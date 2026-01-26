#!/bin/bash
#SBATCH -p gpu41
#SBATCH --gres=gpu:1
#SBATCH -N 1
#SBATCH -J flowpacker

echo "========== Job started at: $(date) =========="
start_time=$(date +%s)

/lustre/grp/cmclab/share/wuxy/anaconda/envs/flowpacker/bin/python demo_scripts/flowpacker_af3score/sampler_pdb_pipe.py $1  --save_dir $2 --use_gt_masks True --csv_file $3



end_time=$(date +%s)
elapsed=$((end_time - start_time))

echo "========== Job finished at: $(date) =========="
echo "========== Total runtime: ${elapsed} seconds =========="
