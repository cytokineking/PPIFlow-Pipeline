#!/bin/bash
#SBATCH -p cpu1,cpu2,cpu_short,fat,hygon,dcu
#SBATCH -n 1
#SBATCH -o per-residue-%j.out
#SBATCH -J per-residue

xml=$1
pdb=$2
log=$3


/lustre/grp/cmclab/share/yuql/tool/rosetta.binary.linux.release-371/main/source/bin/rosetta_scripts.static.linuxgccrelease \
  -parser:protocol $xml \
  -s $pdb \
  -overwrite \
  -ignore_zero_occupancy false > $log
