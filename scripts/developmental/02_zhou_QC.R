# load required packages
library(here)
library(tidyverse)
library(Seurat)
library(scDblFinder)
library(BiocParallel)
library(HGNChelper)
library(Matrix)

# set working directory
setwd(here("data"))

# get data location
dirs <- list.dirs(path = "zhou_2022/",
                  recursive = FALSE, full.names = FALSE)

# import matrix files and create seurat object for each folder in dirs
for(x in dirs) {
  name <- gsub("_matrix", "", x)
  
  cts <- ReadMtx(mtx = paste0("./zhou_2022/", x, "/matrix.mtx.gz"),
                 features = paste0("./zhou_2022/", x, "/features.tsv.gz"),
                 cells = paste0("./zhou_2022/", x, "/barcodes.tsv.gz"))
  
  assign(name, CreateSeuratObject(counts = cts))
}

# merge seurat objects from the same biological replicate together
GW7_sample1_whole_male <- merge(GW7_sample1_lane1_whole_male, GW7_sample1_lane2_whole_male,
                                add.cell.ids = c("lane1", "lane2"), merge.data = FALSE)

GW8_sample1_whole_female <- merge(GW8_sample1_lane1_whole_female, GW8_sample1_lane2_whole_female,
                                  add.cell.ids = c("lane1", "lane2"), merge.data = FALSE)

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
for (name in names(zhou_objects)) {
  info <- str_match(name, "^(GW\\d+)_(sample\\d+)_(whole|anterior|medial|posterior)_(male|female)$")
  
  obj <- zhou_objects[[name]]
  
  obj$gestational.week <- rep(info[2], ncol(obj)) 
  obj$sample.id <- rep(info[3], ncol(obj)) 
  obj$hypothalamus.location <- rep(info[4], ncol(obj)) 
  obj$sex <- rep(info[5], ncol(obj))
  obj$sample <- rep(name, ncol(obj))
  obj$donor <- rep(paste(info[2], info[3], sep = "_"), ncol(obj)) # biological donor, irrespective of anatomical dissection
  obj$barcode <- colnames(obj)
  
  zhou_objects[[name]] <- obj
}

# harmonise gene symbols of each sample objects individually
zhou_symbol_map <- NULL

for (name in names(zhou_objects)) {
  
  message("Harmonising gene symbols for: ", name)
  
  harmonised_result <- harmonise_gene_symbols(zhou_objects[[name]], dataset_name = "zhou")
  
  # replace the original sample object with its corrected version
  zhou_objects[[name]] <- harmonised_result$object
  
  # retain the mapping table
  if (is.null(zhou_symbol_map)) { 
    zhou_symbol_map <- harmonised_result$symbol_map
    } else {
    new_symbols <- !harmonised_result$symbol_map$original_symbol %in%
      zhou_symbol_map$original_symbol
    zhou_symbol_map <- rbind(zhou_symbol_map, harmonised_result$symbol_map[new_symbols, ])
  }
  
  # remove temporary result before processing the next sample
  rm(harmonised_result)
  gc()
}

# save the mapping table for reproducibility
write.csv(zhou_symbol_map, here("csv", "zhou_gene_symbol_mapping.csv"), row.names = FALSE)
table(zhou_symbol_map$changed)

# validate every corrected sample
vapply(zhou_objects,
  function(object) {
    anyDuplicated(rownames(object))
  }, integer(1))

# merged object 
seurat_zhou <- merge(x = zhou_objects[[1]], y = zhou_objects[-1],
                     add.cell.ids = names(zhou_objects), project = "zhou_2022", merge.data = FALSE)

# check metadata
head(seurat_zhou@meta.data)

# check if all sample info are in the metadata columns
unique(seurat_zhou@meta.data$gestational.week) 
unique(seurat_zhou@meta.data$sample) 
unique(seurat_zhou@meta.data$hypothalamus.location)
unique(seurat_zhou@meta.data$sex)

# check cell number before QC
seurat_zhou # 295,985 cells

# calculate mitochondiral RNA percentage
seurat_zhou$percent.mt <- PercentageFeatureSet(seurat_zhou, pattern = "^MT-")

# calculate ribosomal RNA percentage
seurat_zhou$percent.rb <- PercentageFeatureSet(seurat_zhou, pattern = "^RP[SL]")

# calculate haemoglobin RNA percentage
seurat_zhou$percent.hb <- PercentageFeatureSet(seurat_zhou, features = hb_genes)

# check QC columns have been added to metadata
head(seurat_zhou@meta.data)

# filter out cells expressing more than 200 and less than 4000 genes, 10% mtDNA and 5% hbDNA
seurat_zhou <- subset(seurat_zhou,
                      subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & 
                        nCount_RNA > 1000 & nCount_RNA < 15000 &
                        percent.mt < 10 & percent.hb < 5)

# check cell number after initial QC
seurat_zhou # 230,407 cells

# create SingleCellExperiment from seurat object for doublet removal
seurat_zhou <- JoinLayers(seurat_zhou)
sce_zhou <- as.SingleCellExperiment(seurat_zhou)

# doublet removal using scDblFinder
sce_zhou <- scDblFinder(sce_zhou, samples = sce_zhou$sample, BPPARAM = bp)

# check how many singlets and doublets there are
table(sce_zhou$scDblFinder.class) # 33,875 doublets

# convert sce back into seurat object
doublets_zhou <- as.Seurat(sce_zhou, counts = "counts", data = NULL)

# remove doublets
seurat_zhou <- subset(doublets_zhou, subset = scDblFinder.class  == "singlet")

# check how many cells after doublet removal
seurat_zhou # 196,532 cells

# save clean zhou seurat obj
saveRDS(seurat_zhou, here("r_obj", "seurat_zhou.rds"))