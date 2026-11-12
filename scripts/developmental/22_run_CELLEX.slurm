#!/bin/bash -l

#SBATCH --job-name=cellex
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/cellex_%j.out
#SBATCH --error=logs/cellex_%j.err

set -euo pipefail

# activate CELLECT conda environment
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate cellect

# run CELLEX expression specificity analysis
python run_cellex.py