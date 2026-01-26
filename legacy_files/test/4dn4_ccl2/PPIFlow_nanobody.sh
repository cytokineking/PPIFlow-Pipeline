#!/bin/bash
#SBATCH -p gpu41,gpu43,gpu53,gpu51
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -N 1
#SBATCH -o logs/%j_%x.out
#SBATCH -J PPIFlow_nanobody
export MPLBACKEND=Agg

echo "Job started at: $(date)"
start_time=$(date +%s)

ANTIGEN_PDB="/lustre/grp/cmclab/qinxy/antibody_ppiflow/demo/input/4dn4.pdb"
FRAMEWORK_PDB="/lustre/grp/cmclab/qinxy/antibody_ppiflow/demo/input/7xl0_nanobody_framework.pdb"
SPECIFIED_HOTSPOTS="M28,M39,M55,M61"
CDR_LENGTH="CDRH1,8-8,CDRH2,8-8,CDRH3,9-21"
SAMPLES_PER_TARGET="10"
MODEL_WEIGHTS="/lustre/grp/cmclab/share/huangxk/nanobody_1103/epoch=124-step=13000.ckpt"
OUTPUT_DIR="/lustre/grp/cmclab/qinxy/antibody_ppiflow/demo/test/4dn4_ccl2"
NAME=$(basename $OUTPUT_DIR)

ANTIGEN_CHAIN="M"
HEAVY_CHAIN="A"
LIGHT_CHAIN="B"
CONFIG="/lustre/grp/cmclab/qinxy/software/ppiflow_nanobody/configs/inference_nanobody.yaml"

/lustre/grp/cmclab/qinxy/new_miniconda3/envs/fm56/bin/python sample_antibody_nanobody.py \
    --antigen_pdb $ANTIGEN_PDB \
    --framework_pdb $FRAMEWORK_PDB \
    --antigen_chain $ANTIGEN_CHAIN \
    --heavy_chain $HEAVY_CHAIN \
    --specified_hotspots $SPECIFIED_HOTSPOTS \
    --cdr_length $CDR_LENGTH \
    --samples_per_target $SAMPLES_PER_TARGET \
    --config $CONFIG \
    --model_weights $MODEL_WEIGHTS \
    --output_dir $OUTPUT_DIR \
    --name $NAME

end_time=$(date +%s)
echo "Job ended at: $(date)"
elapsed=$((end_time - start_time))
echo "Elapsed time: $elapsed seconds"
