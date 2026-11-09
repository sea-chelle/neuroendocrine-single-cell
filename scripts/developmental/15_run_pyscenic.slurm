#!/bin/bash -l

#SBATCH --job-name=pySCENIC
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --gres=tmpfs:50G
#SBATCH --chdir=/home/sejj745/Scratch/pySCENIC
#SBATCH --output=/home/sejj745/Scratch/pySCENIC/logs/%x_%j.out
#SBATCH --error=/home/sejj745/Scratch/pySCENIC/logs/%x_%j.err

set -euo pipefail

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

module load ucl-stack/2026-07
module load default-modules/2026-07

ROOT="/home/sejj745/Scratch/pySCENIC"
INPUT="${ROOT}/input"
RESOURCES="${ROOT}/resources"
OUTPUT="${ROOT}/output"
CONTAINER="${ROOT}/container/aertslab-pyscenic-0.12.1.sif"
CORES="${SLURM_CPUS_PER_TASK}"

mkdir -p "${OUTPUT}"

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Cores: ${CORES}"
echo "Memory allocation: ${SLURM_MEM_PER_NODE:-unknown} MB"

required_files=(
  "${INPUT}/pyscenic_inference.loom"
  "${INPUT}/pyscenic_scoring.loom"
  "${RESOURCES}/allTFs_hg38.txt"
  "${RESOURCES}/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
  "${RESOURCES}/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
  "${RESOURCES}/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
  "${CONTAINER}"
)

for file in "${required_files[@]}"; do
  if [[ ! -s "${file}" ]]; then
    echo "ERROR: missing or empty file: ${file}" >&2
    exit 1
  fi
done

echo "All required files were found."


# stage 1: GRNBoost2

ADJACENCIES="${OUTPUT}/pyscenic_adjacencies.tsv"

if [[ ! -s "${ADJACENCIES}" ]]; then

  echo "Starting GRNBoost2: $(date)"

  apptainer exec \
    --bind "${ROOT}:/data" \
    "${CONTAINER}" \
    pyscenic grn \
      --method grnboost2 \
      --num_workers "${CORES}" \
      --seed 123 \
      --sparse \
      --output /data/output/pyscenic_adjacencies.tsv \
      /data/input/pyscenic_inference.loom \
      /data/resources/allTFs_hg38.txt

  echo "GRNBoost2 completed: $(date)"

else

  echo "Existing adjacency file found; skipping GRNBoost2."

fi

if [[ ! -s "${ADJACENCIES}" ]]; then
  echo "ERROR: GRNBoost2 did not produce an adjacency file." >&2
  exit 1
fi


# stage 2: motif enrichment and pruning

REGULONS="${OUTPUT}/pyscenic_regulons.csv"

if [[ ! -s "${REGULONS}" ]]; then

  echo "Starting motif pruning: $(date)"

  apptainer exec \
    --bind "${ROOT}:/data" \
    "${CONTAINER}" \
    pyscenic ctx \
      /data/output/pyscenic_adjacencies.tsv \
      /data/resources/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
      /data/resources/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
      --annotations_fname /data/resources/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
      --expression_mtx_fname /data/input/pyscenic_inference.loom \
      --mode custom_multiprocessing \
      --num_workers "${CORES}" \
      --output /data/output/pyscenic_regulons.csv

  echo "Motif pruning completed: $(date)"

else

  echo "Existing regulon file found; skipping motif pruning."

fi

if [[ ! -s "${REGULONS}" ]]; then
  echo "ERROR: motif pruning did not produce a regulon file." >&2
  exit 1
fi


# stage 3: AUCell scoring

AUC_LOOM="${OUTPUT}/pyscenic_scoring_auc.loom"

if [[ ! -s "${AUC_LOOM}" ]]; then

  echo "Starting AUCell: $(date)"

  apptainer exec \
    --bind "${ROOT}:/data" \
    "${CONTAINER}" \
    pyscenic aucell \
      --num_workers "${CORES}" \
      --seed 123 \
      --output /data/output/pyscenic_scoring_auc.loom \
      /data/input/pyscenic_scoring.loom \
      /data/output/pyscenic_regulons.csv

  echo "AUCell completed: $(date)"

else

  echo "Existing AUCell loom found; skipping AUCell."

fi

if [[ ! -s "${AUC_LOOM}" ]]; then
  echo "ERROR: AUCell did not produce its output loom." >&2
  exit 1
fi

# stage 4: export AUCell activity matrix

AUC_CSV="${OUTPUT}/pyscenic_auc_matrix.csv.gz"

if [[ ! -s "${AUC_CSV}" ]]; then

  echo "Exporting AUCell activity matrix: $(date)"

  apptainer exec \
    --bind "${ROOT}:/data" \
    "${CONTAINER}" \
    python - <<'PY'
import h5py
import pandas as pd

# open pySCENIC AUCell loom
with h5py.File("/data/output/pyscenic_scoring_auc.loom", "r") as f:

    # extract cell identifiers
    cell_ids = [
        x.decode() if isinstance(x, bytes) else x
        for x in f["col_attrs"]["CellID"][:]
    ]

    # extract regulon activity scores
    auc = f["col_attrs"]["RegulonsAUC"][:]

    # convert to dataframe
    auc_df = pd.DataFrame.from_records(auc)
    auc_df.index = cell_ids
    auc_df.index.name = "cell"

    print("AUCell matrix dimensions:", auc_df.shape)

    # export compressed CSV
    auc_df.to_csv(
        "/data/output/pyscenic_auc_matrix.csv.gz",
        compression="gzip"
    )

print("AUCell CSV export completed.")
PY

  echo "AUCell matrix exported: $(date)"

else

  echo "Existing AUCell CSV found; skipping export."

fi

if [[ ! -s "${AUC_CSV}" ]]; then
  echo "ERROR: AUCell CSV export failed." >&2
  exit 1
fi

echo
echo "pySCENIC completed successfully."
echo "Job finished: $(date)"
echo

ls -lh "${OUTPUT}"