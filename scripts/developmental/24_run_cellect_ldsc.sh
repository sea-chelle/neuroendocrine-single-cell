#!/bin/bash -l
#SBATCH --job-name=cellect_LDSC
#SBATCH --partition=myriad
#SBATCH --qos=highcpu48h
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --gres=tmpfs:50G
#SBATCH --output=/home/sejj745/Scratch/CELLECT_NEN/logs/cellect_ldsc_%j.out
#SBATCH --error=/home/sejj745/Scratch/CELLECT_NEN/logs/cellect_ldsc_%j.err

set -euo pipefail

# activate CELLECT environment
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate cellect

# run CELLECT-LDSC
cd "$HOME/Scratch/CELLECT"

snakemake \
  --use-conda \
  --conda-frontend conda \
  --cores 4 \
  -s cellect-ldsc.snakefile \
  --configfile "$HOME/Scratch/CELLECT_NEN/config_ldsc.yml"
