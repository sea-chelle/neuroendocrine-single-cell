# load required packages
library(here)
library(Seurat)
library(scDblFinder)
library(BiocParallel)
library(HGNChelper)
library(Matrix)
library(dplyr)

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

# create a function that maps symbols and collapses duplicate rows
harmonise_gene_symbols <- function(object, dataset_name, assay = "RNA") {
  
  DefaultAssay(object) <- assay
  
  # join count layers if the object contains more than one
  all_layers <- Layers(object[[assay]])
  count_layers <- grep("^counts", all_layers, value = TRUE)
  
  if (length(count_layers) == 0) {
    stop("No raw counts layer found in the RNA assay.")
  }
  
  if (length(count_layers) > 1) {
    object <- JoinLayers(object, assay = assay, layers = count_layers, new = "counts")
  }
  
  # extract the raw counts
  counts <- LayerData(object, assay = assay, layer = "counts")
  original_symbols <- rownames(counts)
  
  # map old aliases to approved HGNC symbols
  symbol_check <- checkGeneSymbols(original_symbols, species = "human", unmapped.as.na = FALSE, expand.ambiguous = FALSE)
  approved_symbols <- symbol_check$Suggested.Symbol
  
  # safety fallback for any missing mappings
  approved_symbols[is.na(approved_symbols) | approved_symbols == ""] <- original_symbols[is.na(approved_symbols) | 
                                                                                           approved_symbols == ""]
  
  symbol_map <- data.frame(dataset = dataset_name, 
                           original_symbol = original_symbols, 
                           approved_symbol = approved_symbols,
                           changed = original_symbols != approved_symbols,
                           stringsAsFactors = FALSE)
  
  # create a sparse mapping matrix: old gene rows -> approved gene rows
  unique_symbols <- unique(approved_symbols)
  
  collapse_matrix <- sparseMatrix(i = match(approved_symbols, unique_symbols),
                                  j = seq_along(approved_symbols),
                                  x = 1,
                                  dims = c(length(unique_symbols), length(approved_symbols)),
                                  dimnames = list(unique_symbols, original_symbols))
  
  # sum aliases that map to the same approved symbol
  corrected_counts <- collapse_matrix %*% counts
  
  # confirm that total counts per cell have not changed
  stopifnot(all(Matrix::colSums(counts) == Matrix::colSums(corrected_counts)))
  
  # preserve metadata, but recalculate RNA QC variables
  metadata <- object[[]]
  metadata <- metadata[colnames(corrected_counts), , drop = FALSE]
  
  metadata$nCount_RNA <- NULL
  metadata$nFeature_RNA <- NULL
  metadata$percent.mt <- NULL
  
  corrected_object <- CreateSeuratObject(counts = corrected_counts, meta.data = metadata, project = dataset_name)
  
  corrected_object[["percent.mt"]] <- PercentageFeatureSet(corrected_object, pattern = "^MT-")
  
  return(list(object = corrected_object, symbol_map = symbol_map))
}

# harmonise gene symbols
herb_harmonised <- harmonise_gene_symbols(seurat_herb, dataset_name = "herb")

# retain the corrected Seurat object
seurat_herb <- herb_harmonised$object

# save the mapping table for reproducibility
write.csv(herb_harmonised$symbol_map, here("csv", "herb_gene_symbol_mapping.csv"), row.names = FALSE)

# validate harmonisation
table(herb_harmonised$symbol_map$changed)
anyDuplicated(rownames(seurat_herb))
Layers(seurat_herb[["RNA"]])

# create gestational-week mapping
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

# add gestational week
seurat_herb$gestational.week <- unname(gw_mapping[as.character(seurat_herb$sample)])

# check cell number before QC
seurat_herb # 51,554 cells

# calculate mitochondrial DNA percentage
seurat_herb$percent.mt <- PercentageFeatureSet(seurat_herb, pattern = "^MT-")

# calculate ribosomal DNA percentage
seurat_herb$percent.rb <- PercentageFeatureSet(seurat_herb, pattern = "^RP[SL]")

# calculate haemoglobin DNA percentage
hb_genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2", "HBE1", "HBM", "HBQ1", "HBZ")
seurat_herb$percent.hb <- PercentageFeatureSet(seurat_herb, features = hb_genes)

# filter out cells expressing more than 10% mtDNA and 5% hbDNA and keep cells between 200-4000 total reads
seurat_herb <- subset(seurat_herb, 
                      subset = nFeature_RNA >= 200 & nFeature_RNA <= 4000 &
                        percent.mt < 10 &
                        percent.hb < 5)

# check how many cells after initial QC
seurat_herb # 44912 cells

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
saveRDS(seurat_herb, here("r_obj", "seurat_herb.rds"))