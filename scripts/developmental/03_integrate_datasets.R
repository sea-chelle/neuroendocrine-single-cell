# load required packages
library(Seurat)
library(tidyverse)
library(SingleCellExperiment)
library(scDblFinder)
library(BiocParallel)
library(harmony)

# set working directory
setwd(here("r_obj"))

# merge herb and zhou datasets
seurat_dev <- merge(seurat_herb, seurat_zhou, add.cell.ids = c("herb", "zhou"), merge.data = FALSE)

# perform standard seurat workflow
seurat_dev <- NormalizeData(seurat_dev)
seurat_dev <- FindVariableFeatures(seurat_dev, selection.method = "vst", nfeatures = 2000)

# plot top 10 most variable genes among cells
LabelPoints(plot = VariableFeaturePlot(seurat_dev),
            points = head(VariableFeatures(seurat_dev), 10),
            repel = TRUE, xnudge = 0, ynudge = 0)

# calculate cell cycle scores
seurat_dev <- CellCycleScoring(object = seurat_dev, 
                               g2m.features = cc.genes$g2m.genes,
                               s.features = cc.genes$s.genes)

# check cell cycle scores have been added to metadata
head(seurat_dev@meta.data)

# continue with standard seurat workflow
seurat_dev <- ScaleData(seurat_dev, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
seurat_dev <- RunPCA(seurat_dev, features = VariableFeatures(object = seurat_dev))
ElbowPlot(seurat_dev)

# extract cell names
cell_names <- rownames(seurat_dev@meta.data)

# extract the dataset from the cell name and put in a new column called dataset in metadata
seurat_dev$dataset <- sub("_.*", "", cell_names)

# check that the dataset column has been added
head(seurat_dev@meta.data)
table(seurat_dev$dataset)
  
# run harmony for data integration and batch correction
seurat_dev <- RunHarmony(seurat_dev, group.by.vars = c("dataset", "sample"))

# create donor column
seurat_dev$donor <- as.character(seurat_dev$sample)

# collapse Zhou dissections from the same biological donor
seurat_dev$donor[seurat_dev$dataset == "zhou"] <- sub("_(anterior|medial|posterior|whole)_(male|female)$",
                                                      "",
                                                      seurat_dev$sample[seurat_dev$dataset == "zhou"])

# save integrated seurat object
saveRDS(seurat_dev, "seurat_dev.rds")

# remove unecessary objects
rm("seurat_herb", "seurat_zhou", "cell_names")