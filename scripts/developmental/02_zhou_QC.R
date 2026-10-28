# load required packages
library(here)
library(tidyverse)
library(Seurat)
library(scDblFinder)

# set working directory
setwd(here("data"))

# get data location
dirs <- list.dirs(path = "zhou_2022/",
                  recursive = FALSE, full.names = FALSE)

# importing matrix files and creating seurat object for each folder in dirs
for(x in dirs) {
  name <- gsub("_matrix", "", x)
  
  cts <- ReadMtx(mtx = paste0("./zhou_2022/", x, "/matrix.mtx.gz"),
                 features = paste0("./zhou_2022/", x, "/features.tsv.gz"),
                 cells = paste0("./zhou_2022/", x, "/barcodes.tsv.gz"))
  
  assign(name, CreateSeuratObject(counts = cts))
}

# merge seurat objects from the same biological replicate together
GW7_sample1_whole_male <- merge(GW7_sample1_lane1_whole_male, GW7_sample1_lane2_whole_male)
GW8_sample1_whole_female <- merge(GW8_sample1_lane1_whole_female, GW8_sample1_lane2_whole_female)
rm(list = c("GW7_sample1_lane1_whole_male", "GW7_sample1_lane2_whole_male", 
            "GW8_sample1_lane1_whole_female", "GW8_sample1_lane2_whole_female"))

# get names of all zhou samples
sample_names <- c("GW10_sample1_whole_male", "GW12_sample1_whole_male", "GW12_sample2_whole_male", 
                  "GW15_sample1_anterior_female", "GW15_sample1_medial_female", "GW15_sample1_posterior_female", 
                  "GW18_sample1_anterior_male", "GW18_sample1_medial_male", "GW18_sample1_posterior_male", 
                  "GW18_sample2_anterior_male", "GW18_sample2_medial_male", "GW18_sample2_posterior_male",
                  "GW20_sample1_anterior_male", "GW20_sample1_medial_male", "GW20_sample1_posterior_male", 
                  "GW7_sample1_whole_male", "GW8_sample1_whole_female")

# collect the original sample objects in the order listed
zhou_objects <- mget(sample_names)

# add metadata from each object's known name
for (nm in names(zhou_objects)) {
  info <- str_match(nm, "^(GW\\d+)_(sample\\d+)_(whole|anterior|medial|posterior)_(male|female)$")
  
  obj <- zhou_objects[[nm]]
  obj$gestational.week <- rep(info[2], ncol(obj))
  obj$sample.id <- rep(info[3], ncol(obj))
  obj$hypothalamus.location <- rep(info[4], ncol(obj))
  obj$sex <- rep(info[5], ncol(obj))
  obj$sample <- rep(nm, ncol(obj))
  obj$barcode <- colnames(obj)
  zhou_objects[[nm]] <- obj
}

# merged object 
seurat_zhou <- merge(x = zhou_objects[[1]], y = zhou_objects[-1],
                     add.cell.ids = names(zhou_objects), project = "zhou_2022")

# check metadata
head(seurat_zhou@meta.data)

# check if all sample info are in the metadata columns
unique(seurat_zhou@meta.data$gestational.week) # checking gestational.week column
unique(seurat_zhou@meta.data$sample) #checking sample column
unique(seurat_zhou@meta.data$hypothalamus.location) #checking hypothalamus.location column
unique(seurat_zhou@meta.data$sex) # checking sex column

# check how many cells before QC
seurat_zhou # 295,985 cells

# calculate mitochondiral DNA percentage
seurat_zhou$percent.mt <- PercentageFeatureSet(seurat_zhou, pattern = "^MT-")

# calculate ribosomal DNA percentage
seurat_zhou$percent.rb <- PercentageFeatureSet(seurat_zhou, pattern = "^RP[SL]")

# calculate haemoglobin DNA percentage
seurat_zhou$percent.hb <- PercentageFeatureSet(seurat_zhou, features = hb_genes)

# check QC columns have been added to metadata
head(seurat_zhou@meta.data)

# filter out cells expressing more than 200 and less than 4000 genes, 10% mtDNA and 5% hbDNA
seurat_zhou <- subset(seurat_zhou,
                      subset = nFeature_RNA > 200 &
                        nFeature_RNA < 4000 &
                        nCount_RNA > 1000 &
                        nCount_RNA < 15000 &
                        percent.mt < 10 &
                        percent.hb < 5)

# check how many cells after initial QC
seurat_zhou # 230,382 cells

# create SingleCellExperiment from seurat object for doublet removal
seurat_zhou <- JoinLayers(seurat_zhou) # needed for merged seurat objects only
sce_zhou <- as.SingleCellExperiment(seurat_zhou)

# doublet removal using scDblFinder
sce_zhou <- scDblFinder(sce_zhou, samples = sce_zhou$sample)

# check how many singlets and doublets there are
table(sce_zhou$scDblFinder.class) # 13,567 doublets

# convert sce back into seurat object
doublets_zhou <- as.Seurat(sce_zhou, counts = "counts", data = NULL)

# remove doublets
seurat_zhou <- subset(doublets_zhou, subset = scDblFinder.class  == "singlet")

# check how many cells after doublet removal
seurat_zhou # 213,851 cells

# set workding directory to plots
setwd("/Volumes/ritd-ag-project-rd01lx-nhhan20/Michelle/analysed_data/neuroendocrine/plots")

# plot qc violin plots
VlnPlot(seurat_zhou, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb"), ncol = 2, pt.size = 0)
# ggsave("zhou_qc.png", height = 10, width = 5)

# save seurat object
setwd("/Volumes/ritd-ag-project-rd01lx-nhhan20/Michelle/analysed_data/neuroendocrine/data/zhou_2022/")
saveRDS(seurat_zhou, "seurat_zhou.rds")
