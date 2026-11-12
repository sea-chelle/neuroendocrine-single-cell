# load required packages
library(here)
library(readxl)
library(tidyverse)
library(Seurat)

# set working directory
setwd(here("data", "tadross_2025"))

# import tadross 2025 seurat object (inc. 3 samples from herb/siletti dataset)
seurat_adult <- readRDS("human_HYPOMAP_snRNASeq.rds")

# load supplementary tables
excel_sheets("41586_2024_8504_MOESM3_ESM.xlsx")

# import cluster annotations and spatial assignments
adult_annotations <- read_excel("41586_2024_8504_MOESM3_ESM.xlsx", sheet = "6_cluster_tree_annotations")
adult_spatial <- read_excel("41586_2024_8504_MOESM3_ESM.xlsx", sheet = "12_region_abundance_assignment")

# inspect columns
colnames(adult_annotations)
colnames(adult_spatial)

# extract C4 cluster annotations
adult_C4_annotations <- adult_annotations %>%
  filter(clusterlevel == "C4") %>%
  select(cluster_id, cluster_name, top_markers, ncells, Female, Male)

# identify candidate NEN clusters using exact marker names
adult_C4_NEN <- adult_C4_annotations %>%
  filter(str_detect(cluster_name, regex("\\b(AVP|OXT|CRH|TRH|GHRH|GNRH1|SST|TH|DDC)\\b")))

adult_C4_NEN %>%
  select(cluster_name, top_markers, ncells, Female, Male)  %>%
  print(n = Inf)

# select PVN neuropeptide-associated clusters
adult_PVN_clusters <- adult_C4_NEN %>%
  filter(str_detect(cluster_name, "\\b(AVP|OXT|CRH|TRH)\\b")) %>%
  pull(cluster_name)

# PVN neuroendocrine co-markers
PVN_markers <- c("AVP", "OXT", "CRH", "TRH", "NR3C1", "SCGN", "AGTR1", "CALB1", "KCNMB4", "FOXP1", "TAC1", "TH", 
                 "NFIB", "THRB", "DDC")

# order PVN clusters by neuropeptide identity
adult_PVN_order <- adult_PVN_clusters[order(match(case_when(grepl("\\bAVP\\b", adult_PVN_clusters) ~ "AVP",
                                                            grepl("\\bOXT\\b", adult_PVN_clusters) ~ "OXT",
                                                            grepl("\\bCRH\\b", adult_PVN_clusters) ~ "CRH",
                                                            grepl("\\bTRH\\b", adult_PVN_clusters) ~ "TRH"),
                                                  c("AVP", "OXT", "CRH", "TRH")))]

# plot clusters in the specified order
DotPlot(subset(seurat_adult, subset = C4_named %in% adult_PVN_clusters), features = PVN_markers, group.by = "C4_named",  assay = "RNA") +
  scale_y_discrete(limits = rev(adult_PVN_order)) +
  RotatedAxis()
