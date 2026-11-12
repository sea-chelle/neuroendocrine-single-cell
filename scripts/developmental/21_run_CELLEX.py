# import required packages
import os
import pandas as pd
import cellex
from scipy.io import mmread

# define input and output directories
input_dir = "/home/sejj745/Scratch/CELLECT_NEN/input"
output_dir = "/home/sejj745/Scratch/CELLECT_NEN/cellex_output"

# create output directory
os.makedirs(output_dir, exist_ok=True)

# import cell metadata
print("Reading metadata...", flush=True)
metadata = pd.read_csv(f"{input_dir}/metadata.tsv", sep="\t")

# import gene symbols
print("Reading genes...", flush=True)
genes = pd.read_csv(f"{input_dir}/genes.tsv", sep="\t", header=None)[0].astype(str)

# import raw count matrix as sparse matrix
print("Reading count matrix...", flush=True)
counts = mmread(f"{input_dir}/counts.mtx").tocsr()

# check matrix dimensions match genes and cells
print("Matrix shape:", counts.shape, flush=True)
assert counts.shape == (len(genes), len(metadata))

# convert sparse counts to dense float32 dataframe
print("Converting counts to dense float32 CELLEX DataFrame...", flush=True)
expression = pd.DataFrame(counts.astype("float32").toarray(),
                          index=genes,
                          columns=metadata["cell"])

# remove sparse matrix to free memory
del counts

# assign neuronal population annotations to cells
annotation = pd.Series(metadata["cell_type"].values, index=metadata["cell"].values)

# initialise CELLEX expression specificity object
print("Creating CELLEX ESObject...", flush=True)
eso = cellex.ESObject(data=expression,
                      annotation=annotation,
                      dtype="float32",
                      verbose=True)

# calculate cell-type expression specificity
print("Computing expression specificity...", flush=True)
eso.compute(verbose=True)

# export mean expression specificity scores
print("Saving ESmu...", flush=True)
eso.results["esmu"].to_csv(f"{output_dir}/esmu.csv.gz", compression="gzip")

print("CELLEX COMPLETE", flush=True)