# load required packages
library(here)
library(Seurat)
library(tidyverse)
library(Matrix)
library(AnnotationDbi)
library(org.Hs.eg.db)

# set working directory
setwd(here())

# extract raw counts from PVN and ARC neuronal populations
cellex_counts <- GetAssayData(genetics_dev, assay = "RNA", layer = "counts")

# convert gene symbols to Ensembl IDs
cellex_ids <- mapIds(
  org.Hs.eg.db,
  keys = rownames(cellex_counts),
  column = "ENSEMBL",
  keytype = "SYMBOL",
  multiVals = "asNA"
)

# remove unmapped and duplicated gene IDs
keep <- !is.na(cellex_ids) & !duplicated(cellex_ids)
cellex_counts <- cellex_counts[keep, , drop = FALSE]
rownames(cellex_counts) <- unname(cellex_ids[keep])

# check gene mapping
cat("Genes retained:", nrow(cellex_counts), "\n")
stopifnot(!anyDuplicated(rownames(cellex_counts)))

# create metadata for cells in the count matrix
cellex_metadata <- genetics_dev[[]][colnames(cellex_counts),
                                    c("genetics_group", "donor", "sample", "dataset")] %>%
  rownames_to_column("cell") %>%
  rename(cell_type = genetics_group)

# replace hyphens and spaces with underscores
cellex_metadata$cell_type <- gsub("[ -]", "_", as.character(cellex_metadata$cell_type))

# check cell types and alignment
table(cellex_metadata$cell_type, useNA = "ifany")
stopifnot(identical(colnames(cellex_counts), cellex_metadata$cell))

# create directory for CELLEX input files
dir.create("results/CELLECT/input", recursive = TRUE, showWarnings = FALSE)

# export raw counts as sparse matrix
writeMM(cellex_counts, file = "results/CELLECT/input/counts.mtx")

# export genes
write.table(rownames(cellex_counts),
            file = "results/CELLECT/input/genes.tsv",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# export cell metadata
write.table(cellex_metadata,
            file = "results/CELLECT/input/metadata.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)