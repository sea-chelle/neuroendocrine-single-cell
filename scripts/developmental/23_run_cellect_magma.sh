#!/bin/bash -l
#SBATCH --job-name=cellect_MAGMA
#SBATCH --partition=myriad
#SBATCH --qos=highcpu48h
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=/home/sejj745/Scratch/CELLECT_NEN/logs/cellect_magma_%j.out
#SBATCH --error=/home/sejj745/Scratch/CELLECT_NEN/logs/cellect_magma_%j.err

set -euo pipefail

# activate CELLECT environment
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate cellect

# run CELLECT-MAGMA
cd "$HOME/Scratch/CELLECT"

snakemake \
  --use-conda \
  --conda-frontend conda \
  --cores 2 \
  -s cellect-magma.snakefile \
  --configfile "$HOME/Scratch/CELLECT_NEN/config_magma.yml"
