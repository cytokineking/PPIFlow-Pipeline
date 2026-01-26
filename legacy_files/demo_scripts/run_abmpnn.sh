#!/bin/bash
#SBATCH -p gpu41,gpu43,gpu51,gpu53
#SBATCH -N 1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=abmpnn
#SBATCH --output=logs/%x-%j.out

# 设置环境
source /lustre/grp/cmclab/share/anaconda3/etc/profile.d/conda.sh
conda activate mlfold

# 运行参数
ABMPNN_CHECKPOINT="/lustre/grp/cmclab/share/wanghz/Programs/ProteinMPNN-qin/model_weights/"
FOLDER_WITH_PDBS=$1
OUTPUT_DIR=$2
POSITION_FIXED=$3
CHAINS_TO_DESIGN="A"
NUM_SEQ_PER_TARGET=8
SAMPLING_TEMP=0.5

python /lustre/grp/cmclab/share/wanghz/Programs/ProteinMPNN-qin/protein_mpnn_run.py \
        --path_to_model_weights "$ABMPNN_CHECKPOINT" \
        --model_name "abmpnn" \
        --folder_with_pdbs "$FOLDER_WITH_PDBS" \
        --out_folder "$OUTPUT_DIR" \
        --chain_list "$CHAINS_TO_DESIGN" \
        --position_list "$POSITION_FIXED" \
        --num_seq_per_target "$NUM_SEQ_PER_TARGET" \
        --sampling_temp "$SAMPLING_TEMP" \
        --seed 37 \
        --batch_size "$NUM_SEQ_PER_TARGET" \
        --omit_AAs C
