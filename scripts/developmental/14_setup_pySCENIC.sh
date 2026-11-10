# create directories for running pySCENIC on HPC
mkdir -p /pySCENIC/{input,resources,output,logs,scripts,container}

# transfer input files for pySCENIC
rsync -avh --progress \
  pyscenic_inference.loom \
  pyscenic_scoring.loom \
  pyscenic_inference_metadata.csv \
  pyscenic_scoring_metadata.csv \
  pyscenic_genes.txt \
  pyscenic_gene_filter.csv \
  pyscenic_gene_check.csv \
  pyscenic_filter_parameters.csv \
  pySCENIC/input/

# go to resources folder
  cd pySCENIC/resources

# download the human TF list
wget https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt

# download the v10 human motif annotations
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl

# download bothhg38 motif-ranking databases
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather

# pull the pySCENIC container
module load apptainer/1.2.4-1

cd pySCENIC/container

apptainer pull \
  aertslab-pyscenic-0.12.1.sif \
  docker://aertslab/pyscenic:0.12.1

# verify the container
apptainer exec \
  pySCENIC/container/aertslab-pyscenic-0.12.1.sif \
  pyscenic -h
