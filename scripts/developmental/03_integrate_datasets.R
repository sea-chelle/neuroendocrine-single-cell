# load required packages
library(here)
library(tidyverse)
library(Seurat)
library(harmony)

# set working directory
setwd(here("r_obj"))

# add dataset metadata
seurat_herb$dataset <- "herb"
seurat_zhou$dataset <- "zhou"

# merge herb and zhou datasets
seurat_dev <- merge(seurat_herb, seurat_zhou, 
                    add.cell.ids = c("herb", "zhou"), 
                    merge.data = FALSE)

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
  
# run harmony for data integration and batch correction
seurat_dev <- RunHarmony(seurat_dev, group.by.vars = "dataset")

# save integrated seurat object
saveRDS(seurat_dev, "seurat_dev.rds")

# remove everything except seurat_dev
rm(list = setdiff(ls(envir = .GlobalEnv, all.names = TRUE), "seurat_dev"), envir = .GlobalEnv)

# request release of unused memory
gc(full = TRUE)