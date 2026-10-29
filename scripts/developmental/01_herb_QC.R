# load required packages
library(here)
library(Seurat)
library(scDblFinder)
library(BiocParallel)

# set working directory
setwd(here("data"))

# import herb 2023 dataset
load("herb_2022/HypoSamples_PostSCT.rda")

# update seurat objects in HypoNormDat
seurat_herb <- lapply(HypoNormDat, UpdateSeuratObject)

# add sample name as metadata
for (i in seq_along(seurat_herb)) {
  seurat_herb[[i]]$sample <- names(seurat_herb)[i]
}

# merge into one seurat object
seurat_herb <- merge(seurat_herb[[1]], 
                     y = seurat_herb[-1], 
                     add.cell.ids = names(seurat_herb), 
                     project = "herb_2023")

# check all samples are present in seurat object
table(seurat_herb@meta.data$sample)

# extract metadata
metadata_herb <- seurat_herb@meta.data

# create gw mapping vector
gw_mapping <- c("CS13" = "GW6",
                "CS14" = "GW6",
                "CS15" = "GW7",
                "CS22_2_hypo" = "GW8",
                "CS22_hypo" = "GW8",
                "GW16_hypo" = "GW16",
                "GW18_hypo" = "GW18",
                "GW19_hypo" = "GW19",
                "GW20_34_hypo" = "GW20",
                "GW22T_hypo1" = "GW22",
                "GW25_3V_hypo" = "GW25")

# add gestational week column
metadata_herb$gestational.week <- gw_mapping[as.character(metadata_herb$sample)]

# check gestational week column has been added correctly
head(metadata_herb)
unique(metadata_herb$gestational.week)
sum(is.na(metadata_herb$gestational.week)) # check every cell has a value in gestational week column

# add to original seurat metadata
seurat_herb <- AddMetaData(seurat_herb, metadata = metadata_herb)

# check how many cells before QC
seurat_herb # 51,554 cells

# calculate mitochondrial DNA percentage
seurat_herb$percent.mt <- PercentageFeatureSet(seurat_herb, pattern = "^MT-")

# calculate ribosomal DNA percentage
seurat_herb$percent.rb <- PercentageFeatureSet(seurat_herb, pattern = "^RP[SL]")

# calculate haemoglobin DNA percentage
hb_genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2", "HBE1", "HBM", "HBQ1", "HBZ")
seurat_herb$percent.hb <- PercentageFeatureSet(seurat_herb, features = hb_genes)

# check QC columns have been added to metadata
head(seurat_herb@meta.data)

# filter out cells expressing more than 10% mtDNA and 5% hbDNA
seurat_herb <- subset(seurat_herb, subset = percent.mt < 10 & percent.hb < 5)

# check how many cells after initial QC
seurat_herb # 44589 cells

# create SingleCellExperiment from seurat object for doublet removal
sce_herb <- as.SingleCellExperiment(seurat_herb)

# set seed
bp <- SerialParam(RNGseed = 1234)

# doublet removal using scDblFinder
sce_herb <- scDblFinder(sce_herb, samples = sce_herb$sample, BPPARAM = bp)

# check how many singlets and doublets there are
table(sce_herb$scDblFinder.class) # 516 doublets

# convert sce back into seurat object
doublets_herb <- as.Seurat(sce_herb, counts = "counts", data = NULL)

# remove doublets
seurat_herb <- subset(doublets_herb, subset = scDblFinder.class  == "singlet")

# check how many cells after doublet removal
seurat_herb # 44,073 cells

# save clean seurat_herb obj
setwd(here("r_obj"))
saveRDS(seurat_herb, "seurat_herb.rds")