# loading required packages
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(ComplexHeatmap)
library(circlize)

# set working directory to figures
setwd(here("plots"))

# exclude immature neurons for hypothalamic nuclei identification
nuclei_dev <- subset(neuron_dev, subset = cell_type != "immature")

# normalise and scale data as usual
nuclei_dev <- NormalizeData(nuclei_dev)
nuclei_dev <- FindVariableFeatures(nuclei_dev, selection.method = "vst", nfeatures = 2000)

LabelPoints(plot = VariableFeaturePlot(nuclei_dev),
            points = head(VariableFeatures(nuclei_dev), 10),
            repel = TRUE, xnudge = 0, ynudge = 0)

# continue with standard seurat workflow
nuclei_dev <- ScaleData(nuclei_dev, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
nuclei_dev <- RunPCA(nuclei_dev, features = VariableFeatures(object = nuclei_dev))
ElbowPlot(nuclei_dev)

# UMAP dimensionality reduction clusters
nuclei_dev <- FindNeighbors(nuclei_dev, dims = 1:40)
nuclei_dev <- FindClusters(nuclei_dev, resolution = 0.6)
nuclei_dev <- RunUMAP(nuclei_dev, dims = 1:40)

# visualising clusters
DimPlot(nuclei_dev, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
ggsave("UMAP_neuroendocrine_unannotated.png", height = 7, width = 8)

# looking at cluster sizes
table(nuclei_dev@meta.data$seurat_clusters) 

# SCN markers
FeaturePlot(nuclei_dev, features = c("RGS16", "RELN", "ONECUT1", "ONECUT2", "SIX6", "LHX1", "DLX2", "SIX3"))
ggsave("umap_SCN_marker.png", height = 20, width = 30, dpi = 300)

VlnPlot(nuclei_dev, features = c("RGS16", "RELN", "ONECUT1", "ONECUT2", "SIX6", "LHX1", "DLX2", "SIX3"), pt.size = 0)
ggsave("violin_SCN_marker.png", height = 15, width = 35, dpi = 300)

# PVN markers
FeaturePlot(nuclei_dev, features = c("SIM1", "OTP", "POU3F2"))
ggsave("umap_PVN_marker.png", height = 10, width = 20, dpi = 300)

VlnPlot(nuclei_dev, features = c("SIM1", "OTP", "POU3F2"), pt.size = 0)
ggsave("violin_PVN_marker.png", height = 5, width = 30, dpi = 300)

# PO markers
FeaturePlot(nuclei_dev, features = c("ARX", "NFIX", "GAL"))
ggsave("umap_PO_marker.png", height = 10, width = 20, dpi = 300)

VlnPlot(nuclei_dev, features = c("ARX", "NFIX", "GAL"), pt.size = 0)
ggsave("violin_PO_marker.png", height = 5, width = 30, dpi = 300)

# AH markers
FeaturePlot(nuclei_dev, features = c("SIX3", "DLX1", "DLX2", "DLX5", "DLX6", "ACHE"))
ggsave("umap_AH_marker.png", height = 20, width = 15, dpi = 300)

VlnPlot(nuclei_dev, features = c("SIX3", "DLX1", "DLX2", "DLX5", "DLX6", "ACHE"), pt.size = 0)
ggsave("violin_AH_marker.png", height = 10, width = 30, dpi = 300)

# TH markers
FeaturePlot(nuclei_dev, features = c("ISL1", "MEIS2", "SST", "DLK1", "TBX3"))
ggsave("umap_TH_marker.png", height = 20, width = 15, dpi = 300)

VlnPlot(nuclei_dev, features = c("ISL1", "MEIS2", "SST", "DLK1", "TBX3"), pt.size = 0)
ggsave("violin_TH_marker.png", height = 10, width = 30, dpi = 300)

# DMH markers
FeaturePlot(nuclei_dev, features = c("PPP1R17", "GPR50", "FOXP2", "PRLR"))
ggsave("umap_DMH_marker.png", height = 20, width = 20, dpi = 300)

VlnPlot(nuclei_dev, features = c("PPP1R17", "GPR50", "FOXP2", "PRLR"), pt.size = 0)
ggsave("violin_DMH_marker.png", height = 10, width = 30, dpi = 300)

# VMH markers
FeaturePlot(nuclei_dev, features = c("FEZF1", "NR5A1"))
ggsave("umap_VMH_marker.png", height = 10, width = 20, dpi = 300)

VlnPlot(nuclei_dev, features = c("FEZF1", "NR5A1"), pt.size = 0)
ggsave("violin_VMH_marker.png", height = 5, width = 30, dpi = 300)

# ARC markers
FeaturePlot(nuclei_dev, features = c("TBX3", "ISL1", "AGRP", "POMC", "GHRH", "KISS1", "TAC3", "NPY"))
ggsave("umap_ARC_marker.png", height = 20, width = 25, dpi = 300)

VlnPlot(nuclei_dev, features = c("TBX3", "ISL1", "AGRP", "POMC", "GHRH", "KISS1", "TAC3", "NPY"), pt.size = 0)
ggsave("violin_ARC_marker.png", height = 15, width = 30, dpi = 300)

# LH markers
FeaturePlot(nuclei_dev, features = c("LHX9", "PDYN", "HCRT"))
ggsave("umap_LH_marker.png", height = 20, width = 20, dpi = 300)

VlnPlot(nuclei_dev, features = c("LHX9", "PDYN", "HCRT"), pt.size = 0)
ggsave("violin_LH_marker.png", height = 5, width = 30, dpi = 300)

# PH markers
FeaturePlot(nuclei_dev, features = c("CRABP1", "MEIS2", "CRABP2", "IRX1", "IRX2", "IRX4", "NHLH1"))
ggsave("umap_PH_marker.png", height = 20, width = 20, dpi = 300)

VlnPlot(nuclei_dev, features = c("CRABP1", "MEIS2", "CRABP2", "IRX1", "IRX2", "IRX4", "NHLH1"), pt.size = 0)
ggsave("violin_PH_marker.png", height = 15, width = 30, dpi = 300)

# PMN markers
FeaturePlot(nuclei_dev, features = c("ADCYAP1"))
ggsave("umap_PMN_marker.png", height = 5, width = 5, dpi = 300)

VlnPlot(nuclei_dev, features = c("ADCYAP1"), pt.size = 0)
ggsave("violin_PMN_marker.png", height = 5, width = 10, dpi = 300)

# TMN markers
FeaturePlot(nuclei_dev, features = c("HDC", "SLC18A2", "CBLN4", "TBX3"))
ggsave("umap_TMN_marker.png", height = 20, width = 20, dpi = 300)

VlnPlot(nuclei_dev, features = c("HDC", "SLC18A2", "CBLN4", "TBX3"), pt.size = 0)
ggsave("violin_TMN_marker.png", height = 10, width = 30, dpi = 300)

# MN marker visualisation
FeaturePlot(nuclei_dev, features = c("FOXB1", "RPRM", "CALB1", "LHX5", "PITX2"))
ggsave("umap_MN_marker.png", height = 20, width = 15, dpi = 300)

VlnPlot(nuclei_dev, features = c("FOXB1", "RPRM", "CALB1", "LHX5", "PITX2"), pt.size = 0)
ggsave("violin_MN_marker.png", height = 10, width = 30, dpi = 300)

# MN module score visualisation and calculation
nuclei_dev <- AddModuleScore(nuclei_dev, 
                             features = list(MN = c("FOXB1", "RPRM", "CALB1", "LHX5", "PITX2")), 
                             name = "MN")

FeaturePlot(nuclei_dev, features = "MN1", order = TRUE)
ggsave("umap_MN_module.png", height = 5, width = 5, dpi = 300)
VlnPlot(nuclei_dev, features = "MN1", group.by = "seurat_clusters", pt.size = 0)
ggsave("umap_MN_module.png", height = 5, width = 10, dpi = 300)

avg_MN <- aggregate(MN1 ~ seurat_clusters, data = nuclei_dev@meta.data, FUN = mean)
avg_MN <- avg_MN[order(-avg_MN$MN1), ]
head(avg_MN)

# SMN
# SMN marker visualisation
FeaturePlot(nuclei_dev, features = c("LMX1A", "PITX2", "CALB2", "BARHL1", "BARHL2", "NEUROD2"))
ggsave("umap_SMN_marker.png", height = 20, width = 15, dpi = 300)

VlnPlot(nuclei_dev, features = c("LMX1A", "PITX2", "CALB2", "BARHL1", "BARHL2", "NEUROD2"), pt.size = 0)
ggsave("violin_SMN_marker.png", height = 10, width = 30, dpi = 300)

# SMN module score visualisation and calculation
nuclei_dev <- AddModuleScore(nuclei_dev, 
                             features = list(SMN = c("LMX1A", "PITX2", "CALB2", "BARHL1", "BARHL2", "NEUROD2")), 
                             name = "SMN")

FeaturePlot(nuclei_dev, features = "SMN1", order = TRUE)
ggsave("umap_SMN_module.png", height = 5, width = 5, dpi = 300)
VlnPlot(nuclei_dev, features = "SMN1", group.by = "seurat_clusters", pt.size = 0)
ggsave("umap_SMN_module.png", height = 5, width = 10, dpi = 300)

avg_SMN <- aggregate(SMN1 ~ seurat_clusters, data = nuclei_dev@meta.data, FUN = mean)
avg_SMN <- avg_SMN[order(-avg_SMN$SMN1), ]
head(avg_SMN)

# TE
# TE marker visualisation
FeaturePlot(nuclei_dev, features = c("TBR1", "EOMES", "MEIS2"))
ggsave("umap_TE_marker.png", height = 20, width = 20, dpi = 300)

VlnPlot(nuclei_dev, features = c("TBR1", "EOMES", "MEIS2"), pt.size = 0)
ggsave("violin_TE_marker.png", height = 5, width = 30, dpi = 300)

# TE module score visualisation and calculation
nuclei_dev <- AddModuleScore(nuclei_dev, 
                             features = list(TE = c("TBR1", "EOMES", "MEIS2")), 
                             name = "TE")

FeaturePlot(nuclei_dev, features = "TE1", order = TRUE)
ggsave("umap_TE_module.png", height = 5, width = 5, dpi = 300)
VlnPlot(nuclei_dev, features = "TE1", group.by = "seurat_clusters", pt.size = 0)
ggsave("umap_TE_module.png", height = 5, width = 10, dpi = 300)

avg_TE <- aggregate(TE1 ~ seurat_clusters, data = nuclei_dev@meta.data, FUN = mean)
avg_TE <- avg_TE[order(-avg_TE$TE1), ]
head(avg_TE)

# ZI
# ZI marker visualisation
FeaturePlot(nuclei_dev, features = c("MEIS2", "GABRA5"))
ggsave("umap_ZI_marker.png", height = 10, width = 20, dpi = 300)

VlnPlot(nuclei_dev, features = c("MEIS2"), pt.size = 0)
ggsave("violin_ZI_marker.png", height = 5, width = 30, dpi = 300)

# ZI module score visualisation and calculation
nuclei_dev <- AddModuleScore(nuclei_dev, 
                             features = list(ZI = c("MEIS2", "GABRA5")), 
                             name = "ZI")

FeaturePlot(nuclei_dev, features = "ZI1", order = TRUE)
ggsave("umap_ZI_module.png", height = 5, width = 5, dpi = 300)
VlnPlot(nuclei_dev, features = "ZI1", group.by = "seurat_clusters", pt.size = 0)
ggsave("umap_ZI_module.png", height = 5, width = 10, dpi = 300)

avg_ZI <- aggregate(ZI1 ~ seurat_clusters, data = nuclei_dev@meta.data, FUN = mean)
avg_ZI <- avg_ZI[order(-avg_ZI$ZI1), ]
head(avg_ZI)

# Prethalamus
# prethalamus marker visualisation
FeaturePlot(nuclei_dev, features = c("SP9", "ARX"))
ggsave("umap_prethalamus_marker.png", height = 10, width = 20, dpi = 300)

VlnPlot(nuclei_dev, features = c("SP9", "ARX"), pt.size = 0)
ggsave("violin_prethalamus_marker.png", height = 5, width = 30, dpi = 300)

# prethalamus module score visualisation and calculation
nuclei_dev <- AddModuleScore(nuclei_dev, 
                             features = list(prethalamus = c("SP9", "ARX")), 
                             name = "prethalamus")

FeaturePlot(nuclei_dev, features = "prethalamus1", order = TRUE)
ggsave("umap_prethalamus_module.png", height = 5, width = 5, dpi = 300)
VlnPlot(nuclei_dev, features = "prethalamus1", group.by = "seurat_clusters", pt.size = 0)
ggsave("umap_prethalamus_module.png", height = 5, width = 10, dpi = 300)

avg_prethalamus <- aggregate(prethalamus1 ~ seurat_clusters, data = nuclei_dev@meta.data, FUN = mean)
avg_prethalamus <- avg_prethalamus[order(-avg_prethalamus$prethalamus1), ]
head(avg_prethalamus)

# Thalamus
# thalamus marker visualisation
FeaturePlot(nuclei_dev, features = c("GBX2", "TCF7L2", "GATA2", "GATA3"))
ggsave("umap_thalamus_marker.png", height = 20, width = 20, dpi = 300)

VlnPlot(nuclei_dev, features = c("GBX2", "TCF7L2", "GATA2", "GATA3"), pt.size = 0)
ggsave("violin_thalamus_marker.png", height = 10, width = 30, dpi = 300)

# thalamus module score visualisation and calculation
nuclei_dev <- AddModuleScore(nuclei_dev, 
                             features = list(thalamus = c("GBX2", "TCF7L2", "GATA2", "GATA3")), 
                             name = "thalamus")

FeaturePlot(nuclei_dev, features = "thalamus1", order = TRUE)
ggsave("umap_thalamus_module.png", height = 5, width = 5, dpi = 300)
VlnPlot(nuclei_dev, features = "thalamus1", group.by = "seurat_clusters", pt.size = 0)
ggsave("umap_thalamus_module.png", height = 5, width = 10, dpi = 300)

avg_thalamus <- aggregate(thalamus1 ~ seurat_clusters, data = nuclei_dev@meta.data, FUN = mean)
avg_thalamus <- avg_thalamus[order(-avg_thalamus$thalamus1), ]
head(avg_thalamus)

# find differentially expressed genes from each cell type
nuclei_dev_markers <- FindAllMarkers(nuclei_dev, 
                                     only.pos = TRUE, 
                                     min.pct = 0.25,
                                     min.diff.pct = 0.25,
                                     logfc.threshold = 0.25)

# get top 10 markers for each cluster
top10_markers <- function(x){
  nuclei_dev_markers[nuclei_dev_markers$cluster == x, ] %>% head(n = 10)
  }

# create df of top 10 markers for each cluster
top10_markers_neuroendocrine_df <- map_dfr(0:33, top10_markers)

# save as .csv file
setwd(here("csv"))
write.csv(top10_markers_neuroendocrine_df, "nuclei_dev_top10_markers.csv")

# --------------------------------------------------------------------------------
# renaming clusters
# --------------------------------------------------------------------------------

# renaming cluster numbers to cell types
new_cluster_ids_nuclei <- c("TH", "Thalamus", "PVN", "Unknown", "TH", "Unknown", "AH", "Unknown", "TMN", "PO", "PO",
                            "MN", "Unknown", "VMH", "AH", "ARC", "DMH", "PH", "PVN", "AH", "Unknown",
                            "ARC", "LH", "SCN", "PVN", "PH", "SMN", "ARC", "Prethalamus", "PH")

names(new_cluster_ids_nuclei) <- levels(nuclei_dev)
nuclei_dev <- RenameIdents(nuclei_dev, new_cluster_ids_nuclei)

# visualising clusters after cell-type identification
nuclei_dev@meta.data$"cell_type" <- as.factor(nuclei_dev@active.ident)
DimPlot(nuclei_dev, reduction = "umap", group.by="cell_type", label = TRUE) 

ggsave("UMAP_neuroendocrine_annotated.png", height = 7, width = 8)

table(nuclei_dev@meta.data$cell_type) # looking at cluster sizes

# --------------------------------------------------------------------------------
# neuroendocrine cell type cluster
# --------------------------------------------------------------------------------

# UMAP for neuroendocrine markers
FeaturePlot(nuclei_dev, "CRH")
ggsave("UMAP_CRH.png", height = 5, width = 5)

VlnPlot(nuclei_dev, "CRH", pt.size = 0, group.by = "seurat_clusters")
ggsave("vln_CRH.png", height = 5, width = 10)

FeaturePlot(nuclei_dev, "OXT")
ggsave("UMAP_OXT.png", height = 5, width = 5)

VlnPlot(nuclei_dev, "OXT", pt.size = 0, group.by = "seurat_clusters")
ggsave("vln_OXT.png", height = 5, width = 10)

FeaturePlot(nuclei_dev, "AVP")
ggsave("UMAP_AVP.png", height = 5, width = 5)

VlnPlot(nuclei_dev, "AVP", pt.size = 0, group.by = "seurat_clusters")
ggsave("vln_AVP.png", height = 5, width = 10)

FeaturePlot(nuclei_dev, "TRH")
ggsave("UMAP_TRH.png", height = 5, width = 5)

VlnPlot(nuclei_dev, "TRH", pt.size = 0, group.by = "seurat_clusters")
ggsave("vln_TRH.png", height = 5, width = 10)

FeaturePlot(nuclei_dev, "GNRH1")
ggsave("UMAP_GNRH1.png", height = 5, width = 5)

VlnPlot(nuclei_dev, "GNRH1", pt.size = 0, group.by = "seurat_clusters")
ggsave("vln_GNRH1.png", height = 5, width = 10)

FeaturePlot(nuclei_dev, "GHRH")
ggsave("UMAP_GHRH.png", height = 5, width = 5)

VlnPlot(nuclei_dev, "GHRH", pt.size = 0, group.by = "seurat_clusters")
ggsave("vln_GHRH.png", height = 5, width = 10)

FeaturePlot(nuclei_dev, "SST")
ggsave("UMAP_SST.png", height = 5, width = 5)

VlnPlot(nuclei_dev, "SST", pt.size = 0, group.by = "seurat_clusters")
ggsave("vln_SST.png", height = 5, width = 10)

FeaturePlot(nuclei_dev, c("TH", "SLC6A3"))
ggsave("UMAP_dopa.png", height = 10, width = 20)

VlnPlot(nuclei_dev, c("TH", "SLC6A3"), pt.size = 0, group.by = "seurat_clusters")
ggsave("vln_dopa.png", height = 5, width = 30)

# --------------------------------------------------------------------------------
# PVN neuron cell type cluster
# --------------------------------------------------------------------------------

# subsetting neuronal clusters for subclustering
PVN_dev <- subset(nuclei_dev, idents = c("PVN"))

# normalise and scale data as usual
PVN_dev <- NormalizeData(PVN_dev)
PVN_dev <- FindVariableFeatures(PVN_dev, 
                                selection.method = "vst", 
                                nfeatures = 2000)

top10_PVN_dev <- head(VariableFeatures(PVN_dev), 10)
plot_top10_PVN_dev <- VariableFeaturePlot(PVN_dev)

# remove peptides from highly variable genes
peptides <- c("CRH", "OXT", "AVP", "TRH", "SST", "GHRH", "POMC", "NPY", "AGRP", "KISS1", "TAC3", "GNRH1")
hvgs_pvn <- setdiff(VariableFeatures(PVN_dev), peptides)

# save variable genes
ggsave("PVN_dev_top10_variable.png", height = 5, width = 7)

# continue with standard seurat workflow
PVN_dev <- ScaleData(PVN_dev,
                     features = hvgs_pvn,
                     vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
PVN_dev <- RunPCA(PVN_dev, features = hvgs_pvn)
ElbowPlot(PVN_dev)

# UMAP dimensionality reduction clusters
PVN_dev <- FindNeighbors(PVN_dev, dims = 1:20)
PVN_dev <- FindClusters(PVN_dev, resolution = 0.5)
PVN_dev <- RunUMAP(PVN_dev, dims = 1:20)

# visualising clusters
DimPlot(PVN_dev, reduction = "umap", group.by = "seurat_clusters", label = T)

# looking at cluster sizes
table(PVN_dev@meta.data$seurat_clusters) 
ggsave("UMAP_PVN_unannotated.png", height = 5, width = 6)

# find differentially expressed genes from each cell type
PVN_dev_markers <- FindAllMarkers(PVN_dev, 
                                  only.pos = T, 
                                  min.pct = 0.25,
                                  min.diff.pct = 0.25,
                                  logfc.threshold = 0.25)

# get top 10 markers for each cluster
top10_markers <- function(x){
  PVN_dev_markers[PVN_dev_markers$cluster == x, ] %>% head(n = 20)}

# create df of top 10 markers for each cluster
top10_markers_PVN_df <- map_dfr(0:9, top10_markers)
# save as .csv file
setwd("/Volumes/ritd-ag-project-rd01lx-nhhan20/Michelle/analysed_data/neuroendocrine/csv/")
write.csv(top10_markers_PVN_df, "PVN_dev_top10_markers.csv")

# change working directory back to plots
setwd("/Volumes/ritd-ag-project-rd01lx-nhhan20/Michelle/analysed_data/neuroendocrine/plots")

# identifying neuropeptidergic populations in the PVN
FeaturePlot(PVN_dev, features = c("CRH", "AVP", "OXT", "TRH", "SST", "PENK"))
ggsave("UMAP_PVN_neuropeptides.png", height = 20, width = 15, dpi = 300)

VlnPlot(PVN_dev, features = c("CRH", "AVP", "OXT", "TRH", "SST", "PENK"))
ggsave("Vln_PVN_neuropeptides.png", height = 10, width = 30)

# neuroendocrine co-markers
FeaturePlot(PVN_dev, features = c("NR3C1", "SCGN", "AGTR1", "TAC1", "TH", "ADARB2", "CALB1", "KCNMB4", "FOXP1",
                                  "THRB", "NFIB", "CALB2", "SFRP2"))
ggsave("UMAP_PVN_neuroendocrine.png", height = 20, width = 25, dpi = 300)

VlnPlot(PVN_dev, features = c("NR3C1", "SCGN", "PCSK1", "PCSK2", "AGTR1", "TAC1", "TH", "ADARB2", "CALB1", "KCNMB4", "FOXP1",
                              "THRB", "NFIB", "CALB2", "SFRP2"))
ggsave("vln_PVN_neuroendocrine.png", height = 15, width = 25)

# CRHN module score visualisation and calculation
PVN_dev <- AddModuleScore(PVN_dev, 
                          features = list(CRHN = c("CRH", "NR3C1", "AVP", "SCGN", "AGTR1")), 
                          name = "CRHN")

FeaturePlot(PVN_dev, features = "CRHN1", order = TRUE)
ggsave("umap_CRHN_module.png", height = 5, width = 5, dpi = 300)
VlnPlot(PVN_dev, features = "CRHN1", group.by = "seurat_clusters", pt.size = 0)
ggsave("umap_CRHN_module.png", height = 5, width = 10, dpi = 300)

avg_CRHN <- aggregate(CRHN1 ~ seurat_clusters, data = PVN_dev@meta.data, FUN = mean)
avg_CRHN <- avg_CRHN[order(-avg_CRHN$CRHN1), ]
head(avg_CRHN)

# AVP module score visualisation and calculation
PVN_dev <- AddModuleScore(PVN_dev, 
                          features = list(AVP = c("AVP", "TAC1", "TH", "ADARB2")), 
                          name = "AVP")

FeaturePlot(PVN_dev, features = "AVP1", order = TRUE)
ggsave("umap_AVP_module.png", height = 5, width = 5, dpi = 300)
VlnPlot(PVN_dev, features = "AVP1", group.by = "seurat_clusters", pt.size = 0)
ggsave("umap_AVP_module.png", height = 5, width = 10, dpi = 300)

avg_AVP <- aggregate(AVP1 ~ seurat_clusters, data = PVN_dev@meta.data, FUN = mean)
avg_AVP <- avg_AVP[order(-avg_AVP$AVP1), ]
head(avg_AVP)

# OXT module score visualisation and calculation
PVN_dev <- AddModuleScore(PVN_dev, 
                          features = list(OXT = c("OXT", "CALB1", "KCNMB4", "FOXP1", "ADARB2")), 
                          name = "OXT")

FeaturePlot(PVN_dev, features = "OXT1", order = TRUE)
ggsave("umap_OXT_module.png", height = 5, width = 5, dpi = 300)
VlnPlot(PVN_dev, features = "OXT1", group.by = "seurat_clusters", pt.size = 0)
ggsave("umap_OXT_module.png", height = 5, width = 10, dpi = 300)

avg_OXT <- aggregate(OXT1 ~ seurat_clusters, data = PVN_dev@meta.data, FUN = mean)
avg_OXT <- avg_OXT[order(-avg_OXT$OXT1), ]
head(avg_OXT)

# TRH module score visualisation and calculation
PVN_dev <- AddModuleScore(PVN_dev, 
                          features = list(TRH = c("TRH", "AGTR1", "TRHB", "NFIB")), 
                          name = "TRH")

FeaturePlot(PVN_dev, features = "TRH1", order = TRUE)
ggsave("umap_TRH_module.png", height = 5, width = 5, dpi = 300)
VlnPlot(PVN_dev, features = "TRH1", group.by = "seurat_clusters", pt.size = 0)
ggsave("umap_TRH_module.png", height = 5, width = 10, dpi = 300)

avg_TRH <- aggregate(TRH1 ~ seurat_clusters, data = PVN_dev@meta.data, FUN = mean)
avg_TRH <- avg_TRH[order(-avg_TRH$TRH1), ]
head(avg_TRH)

# --------------------------------------------------------------------------------
# PVN neuron cell type cluster
# --------------------------------------------------------------------------------

# subsetting neuronal clusters for subclustering
PVN_dev <- subset(nuclei_dev, idents = c("PVN"))

# normalise and scale data as usual
PVN_dev <- NormalizeData(PVN_dev)
PVN_dev <- FindVariableFeatures(PVN_dev, 
                                selection.method = "vst", 
                                nfeatures = 2000)

top10_PVN_dev <- head(VariableFeatures(PVN_dev), 10)
plot_top10_PVN_dev <- VariableFeaturePlot(PVN_dev)

# remove peptides from highly variable genes
peptides <- c("CRH", "OXT", "AVP", "TRH", "SST", "GHRH", "POMC", "NPY", "AGRP", "KISS1", "TAC3", "GNRH1")
hvgs_pvn <- setdiff(VariableFeatures(PVN_dev), peptides)

# save variable genes
# ggsave("PVN_dev_top10_variable.png", height = 5, width = 7)

# continue with standard seurat workflow
PVN_dev <- ScaleData(PVN_dev,
                     features = hvgs_pvn,
                     vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
PVN_dev <- RunPCA(PVN_dev, features = hvgs_pvn)
ElbowPlot(PVN_dev)

# UMAP dimensionality reduction clusters
PVN_dev <- FindNeighbors(PVN_dev, dims = 1:20)
PVN_dev <- FindClusters(PVN_dev, resolution = 0.5)
PVN_dev <- RunUMAP(PVN_dev, dims = 1:20)

# visualising clusters
DimPlot(PVN_dev, reduction = "umap", group.by = "seurat_clusters", label = T)

# looking at cluster sizes
table(PVN_dev@meta.data$seurat_clusters) 
# ggsave("UMAP_PVN_unannotated.png", height = 5, width = 6)

# find differentially expressed genes from each cell type
PVN_dev_markers <- FindAllMarkers(PVN_dev, 
                                  only.pos = T, 
                                  min.pct = 0.25,
                                  min.diff.pct = 0.25,
                                  logfc.threshold = 0.25)

# get top 10 markers for each cluster
top10_markers <- function(x){
  PVN_dev_markers[PVN_dev_markers$cluster == x, ] %>% head(n = 20)}

# create df of top 10 markers for each cluster
top10_markers_PVN_df <- map_dfr(0:9, top10_markers)
# save as .csv file
setwd("/Volumes/ritd-ag-project-rd01lx-nhhan20/Michelle/analysed_data/neuroendocrine/csv/")
# write.csv(top10_markers_PVN_df, "PVN_dev_top10_markers.csv")

# change working directory back to plots
setwd("/Volumes/ritd-ag-project-rd01lx-nhhan20/Michelle/analysed_data/neuroendocrine/plots")

# identifying neuropeptidergic populations in the PVN
FeaturePlot(PVN_dev, features = c("CRH", "AVP", "OXT", "TRH", "SST", "PENK", "GNRH1"))
# ggsave("UMAP_PVN_neuropeptides.png", height = 25, width = 15, dpi = 300)

VlnPlot(PVN_dev, features = c("CRH", "AVP", "OXT", "TRH", "SST", "PENK", "GNRH1"))
# ggsave("Vln_PVN_neuropeptides.png", height = 15, width = 30)

# neuroendocrine co-markers
FeaturePlot(PVN_dev, features = c("NR3C1", "SCGN", "AGTR1", "TAC1", "TH", "ADARB2", "CALB1", "KCNMB4", "FOXP1",
                                  "THRB", "NFIB", "CALB2", "SFRP2", "OMP", "UCN3", "GHRH"))
# ggsave("UMAP_PVN_neuroendocrine.png", height = 20, width = 25, dpi = 300)

VlnPlot(PVN_dev, features = c("NR3C1", "SCGN", "PCSK1", "PCSK2", "AGTR1", "TAC1", "TH", "ADARB2", "CALB1", "KCNMB4", "FOXP1",
                              "THRB", "NFIB", "CALB2", "SFRP2","OMP", "UCN3", "GHRH"))
# ggsave("vln_PVN_neuroendocrine.png", height = 20, width = 25)

# CRHN module score visualisation and calculation
PVN_dev <- AddModuleScore(PVN_dev, 
                          features = list(CRHN = c("CRH", "NR3C1", "AVP", "SCGN", "AGTR1")), 
                          name = "CRHN")

FeaturePlot(PVN_dev, features = "CRHN1", order = TRUE)
# ggsave("umap_CRHN_module.png", height = 5, width = 5, dpi = 300)
VlnPlot(PVN_dev, features = "CRHN1", group.by = "seurat_clusters", pt.size = 0)
# ggsave("umap_CRHN_module.png", height = 5, width = 10, dpi = 300)

avg_CRHN <- aggregate(CRHN1 ~ seurat_clusters, data = PVN_dev@meta.data, FUN = mean)
avg_CRHN <- avg_CRHN[order(-avg_CRHN$CRHN1), ]
head(avg_CRHN)

# AVP module score visualisation and calculation
PVN_dev <- AddModuleScore(PVN_dev, 
                          features = list(AVP = c("AVP", "TAC1", "TH", "ADARB2")), 
                          name = "AVP")

FeaturePlot(PVN_dev, features = "AVP1", order = TRUE)
# ggsave("umap_AVP_module.png", height = 5, width = 5, dpi = 300)
VlnPlot(PVN_dev, features = "AVP1", group.by = "seurat_clusters", pt.size = 0)
# ggsave("umap_AVP_module.png", height = 5, width = 10, dpi = 300)

avg_AVP <- aggregate(AVP1 ~ seurat_clusters, data = PVN_dev@meta.data, FUN = mean)
avg_AVP <- avg_AVP[order(-avg_AVP$AVP1), ]
head(avg_AVP)

# OXT module score visualisation and calculation
PVN_dev <- AddModuleScore(PVN_dev, 
                          features = list(OXT = c("OXT", "CALB1", "KCNMB4", "FOXP1", "ADARB2")), 
                          name = "OXT")

FeaturePlot(PVN_dev, features = "OXT1", order = TRUE)
# ggsave("umap_OXT_module.png", height = 5, width = 5, dpi = 300)
VlnPlot(PVN_dev, features = "OXT1", group.by = "seurat_clusters", pt.size = 0)
# ggsave("umap_OXT_module.png", height = 5, width = 10, dpi = 300)

avg_OXT <- aggregate(OXT1 ~ seurat_clusters, data = PVN_dev@meta.data, FUN = mean)
avg_OXT <- avg_OXT[order(-avg_OXT$OXT1), ]
head(avg_OXT)

# TRH module score visualisation and calculation
PVN_dev <- AddModuleScore(PVN_dev, 
                          features = list(TRH = c("TRH", "AGTR1", "TRHB", "NFIB", "GHRH", "NFIB", "OMP", "TAC1", "UCN3")), 
                          name = "TRH")

FeaturePlot(PVN_dev, features = "TRH1", order = TRUE)
# ggsave("umap_TRH_module.png", height = 5, width = 5, dpi = 300)
VlnPlot(PVN_dev, features = "TRH1", group.by = "seurat_clusters", pt.size = 0)
# ggsave("umap_TRH_module.png", height = 5, width = 10, dpi = 300)

avg_TRH <- aggregate(TRH1 ~ seurat_clusters, data = PVN_dev@meta.data, FUN = mean)
avg_TRH <- avg_TRH[order(-avg_TRH$TRH1), ]
head(avg_TRH)

# define list of genes for neuroendocrine identification
genes <- c("CRH", "AVP", "OXT", "TRH", "SST", "GHRH", "PENK", "GNRH1", "NR3C1", "SCGN", "PCSK1", "PCSK2", "AGTR1", 
           "TAC1", "TH", "ADARB2", "CALB1", "KCNMB4", "FOXP1","THRB", "NFIB", "CALB2", "SFRP2","OMP", "UCN3",
           "PNOC", "DSCAM", "SP4", "COCH", "CXCL14", "EBF3", "SEMA3C" , "SEMA4D", "SOX4", "ECEL1", "FOXP2")

# initialise list to store results in 
nen_res_PVN <- list()

# loop through each gene
for (gene in genes) {
  
  # pseudobulk expression data for gene per cluster
  avg_expression <- AverageExpression(PVN_dev,
                                      features = gene,
                                      return.seurat = F,
                                      assays = "RNA",
                                      group.by = "seurat_clusters")$RNA
  
  # convert to df and transpose to have clusters as rows
  avg_expression_df <- as.data.frame(t(avg_expression))
  
  # assign gene name as column name
  colnames(avg_expression_df) <- c(gene)
  
  # add cluster column
  avg_expression_df$cluster <- as.numeric(gsub("g", "", rownames(avg_expression_df)))
  
  # melt df for ggplot
  avg_expression_melt <- melt(avg_expression_df, id.vars = "cluster")
  
  # add column names for df
  colnames(avg_expression_melt) <- c("cluster", "gene", "expression")
  
  # calculate the proportion of cells expressing the gene (expression > 0)
  PVN_dev_data <- FetchData(PVN_dev, vars = c(gene, "seurat_clusters"))
  
  proportion_cells <- PVN_dev_data %>%
    group_by(seurat_clusters) %>%
    summarise(percentage = sum(!!sym(gene) > 0) / n() * 100)
  
  # convert "seurat_clusters" to numeric to match with "cluster" in avg_expression_melt
  proportion_cells$seurat_clusters <- as.numeric(as.character(proportion_cells$seurat_clusters))
  
  # calculate mean expression across all clusters
  mean_expr <- mean(avg_expression_melt$expression) 
  
  # calculate log2 fold change for each cluster relative to the global mean expression
  cluster_avg_expr <- avg_expression_melt %>%
    group_by(cluster) %>%
    summarise(avg_expr = mean(expression))
  
  # calculate log2 fold change (LFC) relative to the global mean expression
  cluster_avg_expr$log2fc <- log2(cluster_avg_expr$avg_expr / mean_expr)
  
  # apply all thresholds for filtering:
  neuroendocrine_enriched_clusters <- proportion_cells %>%
    left_join(cluster_avg_expr, by = c("seurat_clusters" = "cluster")) %>%
    filter(
      percentage > 5 &  # percentage of cells expressing > 5%
        log2fc > 0.25 # lfc > 0.25
    )
  
  # store results for  gene in the list
  nen_res_PVN[[gene]] <- neuroendocrine_enriched_clusters
  
  # average expression plot
  ggplot(avg_expression_melt, aes(x = as.factor(cluster), y = expression, fill = as.factor(cluster))) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic() +
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
          legend.position = "none") +
    labs(title = paste("Average expression of", gene),
         x = "Cluster", y = "Average expression") +
    scale_x_discrete(labels = 0:14) +
    scale_y_continuous(expand = c(0,0)) 
  
  ggsave(paste0(gene, "_avg_exp_PVN.png"), height = 4, width = 5)
  
  # percemtage of cells expressing plot
  ggplot(proportion_cells, aes(x = as.factor(seurat_clusters), y = percentage, fill = as.factor(seurat_clusters))) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic() +
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
          legend.position = "none") +
    labs(title = paste("Percentage of Cells Expressing", gene),
         x = "Cluster", y = "Percentage of Cells Expressing") +
    scale_x_discrete(labels = 0:14) +
    scale_y_continuous(expand = c(0, 0)) 
  
  ggsave(paste0(gene, "_per_exp_PVN.png"), height = 4, width = 5)
}

# --------------------------------------------------------------------------------
# PO neuron cell type cluster
# --------------------------------------------------------------------------------

# subsetting neuronal clusters for subclustering
PO_dev <- subset(nuclei_dev, idents = c("PO"))

# normalise and scale data as usual
PO_dev <- NormalizeData(PO_dev)
PO_dev <- FindVariableFeatures(PO_dev, 
                               selection.method = "vst", 
                               nfeatures = 2000)

top10_PO_dev <- head(VariableFeatures(PO_dev), 10)
plot_top10_PO_dev <- VariableFeaturePlot(PO_dev)

# remove peptides from highly variable genes
peptides <- c("CRH", "OXT", "AVP", "TRH", "SST", "GHRH", "POMC", "NPY", "AGRP", "KISS1", "TAC3", "GNRH1")
hvgs_PO <- setdiff(VariableFeatures(PO_dev), peptides)

# save variable genes
# ggsave("PO_dev_top10_variable.png", height = 5, width = 7)

# continue with standard seurat workflow
PO_dev <- ScaleData(PO_dev,
                    features = hvgs_PO,
                    vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
PO_dev <- RunPCA(PO_dev, features = hvgs_PO)
ElbowPlot(PO_dev)

# UMAP dimensionality reduction clusters
PO_dev <- FindNeighbors(PO_dev, dims = 1:20)
PO_dev <- FindClusters(PO_dev, resolution = 0.5)
PO_dev <- RunUMAP(PO_dev, dims = 1:20)

# visualising clusters
DimPlot(PO_dev, reduction = "umap", group.by = "seurat_clusters", label = T)

# looking at cluster sizes
table(PO_dev@meta.data$seurat_clusters) 
# ggsave("UMAP_PO_unannotated.png", height = 5, width = 6)

# find differentially expressed genes from each cell type
PO_dev_markers <- FindAllMarkers(PO_dev, 
                                 only.pos = T, 
                                 min.pct = 0.25,
                                 min.diff.pct = 0.25,
                                 logfc.threshold = 0.25)

# get top 10 markers for each cluster
top10_markers <- function(x){
  PO_dev_markers[PO_dev_markers$cluster == x, ] %>% head(n = 20)}

# create df of top 10 markers for each cluster
top10_markers_PO_df <- map_dfr(0:9, top10_markers)
# save as .csv file
setwd("/Volumes/ritd-ag-project-rd01lx-nhhan20/Michelle/analysed_data/neuroendocrine/csv/")
# write.csv(top10_markers_PO_df, "PO_dev_top10_markers.csv")

# change working directory back to plots
setwd("/Volumes/ritd-ag-project-rd01lx-nhhan20/Michelle/analysed_data/neuroendocrine/plots")

# identifying neuropeptidergic populations in the PO
FeaturePlot(PO_dev, features = c("CRH", "AVP", "OXT", "TRH", "SST", "PENK", "GNRH1"))
# ggsave("UMAP_PO_neuropeptides.png", height = 25, width = 15, dpi = 300)

VlnPlot(PO_dev, features = c("CRH", "AVP", "OXT", "TRH", "SST", "PENK", "GNRH1"))
# ggsave("Vln_PO_neuropeptides.png", height = 15, width = 30)

# neuroendocrine co-markers
FeaturePlot(PO_dev, features = c("NR3C1", "SCGN", "AGTR1", "TAC1", "TH", "ADARB2", "CALB1", "KCNMB4", "FOXP1",
                                 "THRB", "NFIB", "CALB2", "SFRP2", "OMP", "UCN3", "GHRH"))
# ggsave("UMAP_PO_neuroendocrine.png", height = 20, width = 25, dpi = 300)

VlnPlot(PO_dev, features = c("NR3C1", "SCGN", "PCSK1", "PCSK2", "AGTR1", "TAC1", "TH", "ADARB2", "CALB1", "KCNMB4", "FOXP1",
                             "THRB", "NFIB", "CALB2", "SFRP2","OMP", "UCN3", "GHRH"))
# ggsave("vln_PO_neuroendocrine.png", height = 20, width = 25)

# initialise list to store results in 
nen_res_PO <- list()

# loop through each gene
for (gene in genes) {
  
  # pseudobulk expression data for gene per cluster
  avg_expression <- AverageExpression(PO_dev,
                                      features = gene,
                                      return.seurat = F,
                                      assays = "RNA",
                                      group.by = "seurat_clusters")$RNA
  
  # convert to df and transpose to have clusters as rows
  avg_expression_df <- as.data.frame(t(avg_expression))
  
  # assign gene name as column name
  colnames(avg_expression_df) <- c(gene)
  
  # add cluster column
  avg_expression_df$cluster <- as.numeric(gsub("g", "", rownames(avg_expression_df)))
  
  # melt df for ggplot
  avg_expression_melt <- melt(avg_expression_df, id.vars = "cluster")
  
  # add column names for df
  colnames(avg_expression_melt) <- c("cluster", "gene", "expression")
  
  # calculate the proportion of cells expressing the gene (expression > 0)
  PO_dev_data <- FetchData(PO_dev, vars = c(gene, "seurat_clusters"))
  
  proportion_cells <- PO_dev_data %>%
    group_by(seurat_clusters) %>%
    summarise(percentage = sum(!!sym(gene) > 0) / n() * 100)
  
  # convert "seurat_clusters" to numeric to match with "cluster" in avg_expression_melt
  proportion_cells$seurat_clusters <- as.numeric(as.character(proportion_cells$seurat_clusters))
  
  # calculate mean expression across all clusters
  mean_expr <- mean(avg_expression_melt$expression) 
  
  # calculate log2 fold change for each cluster relative to the global mean expression
  cluster_avg_expr <- avg_expression_melt %>%
    group_by(cluster) %>%
    summarise(avg_expr = mean(expression))
  
  # calculate log2 fold change (LFC) relative to the global mean expression
  cluster_avg_expr$log2fc <- log2(cluster_avg_expr$avg_expr / mean_expr)
  
  # apply all thresholds for filtering:
  neuroendocrine_enriched_clusters <- proportion_cells %>%
    left_join(cluster_avg_expr, by = c("seurat_clusters" = "cluster")) %>%
    filter(
      percentage > 5 &  # percentage of cells expressing > 5%
        log2fc > 0.25 # lfc > 0.25
    )
  
  # store results for  gene in the list
  nen_res_PO[[gene]] <- neuroendocrine_enriched_clusters
  
  # average expression plot
  ggplot(avg_expression_melt, aes(x = as.factor(cluster), y = expression, fill = as.factor(cluster))) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic() +
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
          legend.position = "none") +
    labs(title = paste("Average expression of", gene),
         x = "Cluster", y = "Average expression") +
    scale_x_discrete(labels = 0:14) +
    scale_y_continuous(expand = c(0,0)) 
  
  ggsave(paste0(gene, "_avg_exp_PO.png"), height = 4, width = 5)
  
  # percemtage of cells expressing plot
  ggplot(proportion_cells, aes(x = as.factor(seurat_clusters), y = percentage, fill = as.factor(seurat_clusters))) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic() +
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
          legend.position = "none") +
    labs(title = paste("Percentage of Cells Expressing", gene),
         x = "Cluster", y = "Percentage of Cells Expressing") +
    scale_x_discrete(labels = 0:14) +
    scale_y_continuous(expand = c(0, 0)) 
  
  ggsave(paste0(gene, "_per_exp_PO.png"), height = 4, width = 5)
}

# --------------------------------------------------------------------------------
# anterior neuron cell type cluster
# --------------------------------------------------------------------------------

# subsetting neuronal clusters for subclustering
AH_dev <- subset(nuclei_dev, idents = c("AH"))

# normalise and scale data as usual
AH_dev <- NormalizeData(AH_dev)
AH_dev <- FindVariableFeatures(AH_dev, 
                               selection.method = "vst", 
                               nfeatures = 2000)

top10_AH_dev <- head(VariableFeatures(AH_dev), 10)
plot_top10_AH_dev <- VariableFeaturePlot(AH_dev)

# remove peptides from highly variable genes
peptides <- c("CRH", "OXT", "AVP", "TRH", "SST", "GHRH", "AHMC", "NPY", "AGRP", "KISS1", "TAC3", "GNRH1")
hvgs_AH <- setdiff(VariableFeatures(AH_dev), peptides)

# save variable genes
# ggsave("AH_dev_top10_variable.png", height = 5, width = 7)

# continue with standard seurat workflow
AH_dev <- ScaleData(AH_dev,
                    features = hvgs_AH,
                    vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
AH_dev <- RunPCA(AH_dev, features = hvgs_AH)
ElbowPlot(AH_dev)

# UMAP dimensionality reduction clusters
AH_dev <- FindNeighbors(AH_dev, dims = 1:20)
AH_dev <- FindClusters(AH_dev, resolution = 0.5)
AH_dev <- RunUMAP(AH_dev, dims = 1:20)

# visualising clusters
DimPlot(AH_dev, reduction = "umap", group.by = "seurat_clusters", label = T)

# looking at cluster sizes
table(AH_dev@meta.data$seurat_clusters) 
# ggsave("UMAP_AH_unannotated.png", height = 5, width = 6)

# find differentially expressed genes from each cell type
AH_dev_markers <- FindAllMarkers(AH_dev, 
                                 only.AHs = T, 
                                 min.pct = 0.25,
                                 min.diff.pct = 0.25,
                                 logfc.threshold = 0.25)

# get top 10 markers for each cluster
top10_markers <- function(x){
  AH_dev_markers[AH_dev_markers$cluster == x, ] %>% head(n = 20)}

# create df of top 10 markers for each cluster
top10_markers_AH_df <- map_dfr(0:9, top10_markers)
# save as .csv file
setwd("/Volumes/ritd-ag-project-rd01lx-nhhan20/Michelle/analysed_data/neuroendocrine/csv/")
# write.csv(top10_markers_AH_df, "AH_dev_top10_markers.csv")

# change working directory back to plots
setwd("/Volumes/ritd-ag-project-rd01lx-nhhan20/Michelle/analysed_data/neuroendocrine/plots")

# identifying neuropeptidergic AHpulations in the AH
FeaturePlot(AH_dev, features = c("CRH", "AVP", "OXT", "TRH", "SST", "PENK", "GNRH1"))
# ggsave("UMAP_AH_neuropeptides.png", height = 25, width = 15, dpi = 300)

VlnPlot(AH_dev, features = c("CRH", "AVP", "OXT", "TRH", "SST", "PENK", "GNRH1"))
# ggsave("Vln_AH_neuropeptides.png", height = 15, width = 30)

# neuroendocrine co-markers
FeaturePlot(AH_dev, features = c("NR3C1", "SCGN", "AGTR1", "TAC1", "TH", "ADARB2", "CALB1", "KCNMB4", "FOXP1",
                                 "THRB", "NFIB", "CALB2", "SFRP2", "OMP", "UCN3", "GHRH"))
# ggsave("UMAP_AH_neuroendocrine.png", height = 20, width = 25, dpi = 300)

VlnPlot(AH_dev, features = c("NR3C1", "SCGN", "PCSK1", "PCSK2", "AGTR1", "TAC1", "TH", "ADARB2", "CALB1", "KCNMB4", "FOXP1",
                             "THRB", "NFIB", "CALB2", "SFRP2","OMP", "UCN3", "GHRH"))
# ggsave("vln_AH_neuroendocrine.png", height = 20, width = 25)

# initialise list to store results in 
nen_res_AH <- list()

# loop through each gene
for (gene in genes) {
  
  # pseudobulk expression data for gene per cluster
  avg_expression <- AverageExpression(AH_dev,
                                      features = gene,
                                      return.seurat = F,
                                      assays = "RNA",
                                      group.by = "seurat_clusters")$RNA
  
  # convert to df and transpose to have clusters as rows
  avg_expression_df <- as.data.frame(t(avg_expression))
  
  # assign gene name as column name
  colnames(avg_expression_df) <- c(gene)
  
  # add cluster column
  avg_expression_df$cluster <- as.numeric(gsub("g", "", rownames(avg_expression_df)))
  
  # melt df for ggplot
  avg_expression_melt <- melt(avg_expression_df, id.vars = "cluster")
  
  # add column names for df
  colnames(avg_expression_melt) <- c("cluster", "gene", "expression")
  
  # calculate the proportion of cells expressing the gene (expression > 0)
  AH_dev_data <- FetchData(AH_dev, vars = c(gene, "seurat_clusters"))
  
  proportion_cells <- AH_dev_data %>%
    group_by(seurat_clusters) %>%
    summarise(percentage = sum(!!sym(gene) > 0) / n() * 100)
  
  # calculate mean expression across all clusters
  mean_expr <- mean(avg_expression_melt$expression) 
  
  # convert "seurat_clusters" to numeric to match with "cluster" in avg_expression_melt
  proportion_cells$seurat_clusters <- as.numeric(as.character(proportion_cells$seurat_clusters))
  
  # calculate log2 fold change for each cluster relative to the global mean expression
  cluster_avg_expr <- avg_expression_melt %>%
    group_by(cluster) %>%
    summarise(avg_expr = mean(expression))
  
  # calculate log2 fold change (LFC) relative to the global mean expression
  cluster_avg_expr$log2fc <- log2(cluster_avg_expr$avg_expr / mean_expr)
  
  # apply all thresholds for filtering:
  neuroendocrine_enriched_clusters <- proportion_cells %>%
    left_join(cluster_avg_expr, by = c("seurat_clusters" = "cluster")) %>%
    filter(
      percentage > 5 &  # percentage of cells expressing > 5%
        log2fc > 0.25 # lfc > 0.25
    )
  
  # store results for  gene in the list
  nen_res_AH[[gene]] <- neuroendocrine_enriched_clusters
  
  # average expression plot
  ggplot(avg_expression_melt, aes(x = as.factor(cluster), y = expression, fill = as.factor(cluster))) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic() +
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
          legend.position = "none") +
    labs(title = paste("Average expression of", gene),
         x = "Cluster", y = "Average expression") +
    scale_x_discrete(labels = 0:14) +
    scale_y_continuous(expand = c(0,0)) 
  
  ggsave(paste0(gene, "_avg_exp_AH.png"), height = 4, width = 5)
  
  # percemtage of cells expressing plot
  ggplot(proportion_cells, aes(x = as.factor(seurat_clusters), y = percentage, fill = as.factor(seurat_clusters))) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic() +
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
          legend.position = "none") +
    labs(title = paste("Percentage of Cells Expressing", gene),
         x = "Cluster", y = "Percentage of Cells Expressing") +
    scale_x_discrete(labels = 0:14) +
    scale_y_continuous(expand = c(0, 0)) 
  
  ggsave(paste0(gene, "_per_exp_AH.png"), height = 4, width = 5)
}

# --------------------------------------------------------------------------------
# ARC neuron cell type cluster
# --------------------------------------------------------------------------------

# subsetting neuronal clusters for subclustering
ARC_dev <- subset(nuclei_dev, idents = c("ARC"))

# normalise and scale data as usual
ARC_dev <- NormalizeData(ARC_dev)
ARC_dev <- FindVariableFeatures(ARC_dev, 
                                selection.method = "vst", 
                                nfeatures = 2000)

top10_ARC_dev <- head(VariableFeatures(ARC_dev), 10)
plot_top10_ARC_dev <- VariableFeaturePlot(ARC_dev)

# remove peptides from highly variable genes
peptides <- c("CRH", "OXT", "AVP", "TRH", "SST", "GHRH", "POMC", "NPY", "AGRP", "KISS1", "TAC3", "GNRH1")
hvgs_ARC <- setdiff(VariableFeatures(ARC_dev), peptides)

# save variable genes
# ggsave("ARC_dev_top10_variable.png", height = 5, width = 7)

# continue with standard seurat workflow
ARC_dev <- ScaleData(ARC_dev,
                     features = hvgs_ARC,
                     vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
ARC_dev <- RunPCA(ARC_dev, features = hvgs_ARC)
ElbowPlot(ARC_dev)

# UMAP dimensionality reduction clusters
ARC_dev <- FindNeighbors(ARC_dev, dims = 1:20)
ARC_dev <- FindClusters(ARC_dev, resolution = 0.5)
ARC_dev <- RunUMAP(ARC_dev, dims = 1:20)

# visualising clusters
DimPlot(ARC_dev, reduction = "umap", group.by = "seurat_clusters", label = T)

# looking at cluster sizes
table(ARC_dev@meta.data$seurat_clusters) 
# ggsave("UMAP_ARC_unannotated.png", height = 5, width = 6)

# find differentially expressed genes from each cell type
ARC_dev_markers <- FindAllMarkers(ARC_dev, 
                                  only.pos = T, 
                                  min.pct = 0.25,
                                  min.diff.pct = 0.25,
                                  logfc.threshold = 0.25)

# get top 10 markers for each cluster
top10_markers <- function(x){
  ARC_dev_markers[ARC_dev_markers$cluster == x, ] %>% head(n = 20)}

# create df of top 10 markers for each cluster
top10_markers_ARC_df <- map_dfr(0:9, top10_markers)
# save as .csv file
setwd("/Volumes/ritd-ag-project-rd01lx-nhhan20/Michelle/analysed_data/neuroendocrine/csv/")
write.csv(top10_markers_ARC_df, "ARC_dev_top10_markers.csv")

# change working directory back to plots
setwd("/Volumes/ritd-ag-project-rd01lx-nhhan20/Michelle/analysed_data/neuroendocrine/plots")

# identifying neuropeptidergic populations in the ARC
FeaturePlot(ARC_dev, features = c("CRH", "AVP", "OXT", "TRH", "SST", "PENK", "GNRH1"))
# ggsave("UMAP_ARC_neuropeptides.png", height = 25, width = 15, dpi = 300)

VlnPlot(ARC_dev, features = c("CRH", "AVP", "OXT", "TRH", "SST", "PENK", "GNRH1"))
# ggsave("Vln_ARC_neuropeptides.png", height = 15, width = 30)

# neuroendocrine co-markers
FeaturePlot(ARC_dev, features = c("NR3C1", "SCGN", "AGTR1", "TAC1", "TH", "ADARB2", "CALB1", "KCNMB4", "FOXP1",
                                  "THRB", "NFIB", "CALB2", "SFRP2", "OMP", "UCN3", "GHRH"))
# ggsave("UMAP_ARC_neuroendocrine.png", height = 20, width = 25, dpi = 300)

VlnPlot(ARC_dev, features = c("NR3C1", "SCGN", "PCSK1", "PCSK2", "AGTR1", "TAC1", "TH", "ADARB2", "CALB1", "KCNMB4", "FOXP1",
                              "THRB", "NFIB", "CALB2", "SFRP2","OMP", "UCN3", "GHRH"))
# ggsave("vln_ARC_neuroendocrine.png", height = 20, width = 25)

# define list of genes for neuroendocrine identification
genes <- c("CRH", "AVP", "OXT", "TRH", "SST", "GHRH", "PENK", "GNRH1", "NR3C1", "SCGN", "PCSK1", "PCSK2", "AGTR1", 
           "TAC1", "TH", "ADARB2", "CALB1", "KCNMB4", "FOXP1","THRB", "NFIB", "CALB2", "SFRP2","OMP", "UCN3")

# initialise list to store results in 
nen_res_ARC <- list()

# loop through each gene
for (gene in genes) {
  
  # pseudobSulk expression data for gene per cluster
  avg_expression <- AverageExpression(ARC_dev,
                                      features = gene,
                                      return.seurat = F,
                                      assays = "RNA",
                                      group.by = "seurat_clusters")$RNA
  
  # convert to df and transpose to have clusters as rows
  avg_expression_df <- as.data.frame(t(avg_expression))
  
  # assign gene name as column name
  colnames(avg_expression_df) <- c(gene)
  
  # add cluster column
  avg_expression_df$cluster <- as.numeric(gsub("g", "", rownames(avg_expression_df)))
  
  # melt df for ggplot
  avg_expression_melt <- melt(avg_expression_df, id.vars = "cluster")
  
  # add column names for df
  colnames(avg_expression_melt) <- c("cluster", "gene", "expression")
  
  # calculate the proportion of cells expressing the gene (expression > 0)
  ARC_dev_data <- FetchData(ARC_dev, vars = c(gene, "seurat_clusters"))
  
  proportion_cells <- ARC_dev_data %>%
    group_by(seurat_clusters) %>%
    summarise(percentage = sum(!!sym(gene) > 0) / n() * 100)
  
  # convert "seurat_clusters" to numeric to match with "cluster" in avg_expression_melt
  proportion_cells$seurat_clusters <- as.numeric(as.character(proportion_cells$seurat_clusters))
  
  # calculate mean expression across all clusters
  mean_expr <- mean(avg_expression_melt$expression) 
  
  # calculate log2 fold change for each cluster relative to the global mean expression
  cluster_avg_expr <- avg_expression_melt %>%
    group_by(cluster) %>%
    summarise(avg_expr = mean(expression))
  
  # calculate log2 fold change (LFC) relative to the global mean expression
  cluster_avg_expr$log2fc <- log2(cluster_avg_expr$avg_expr / mean_expr)
  
  # apply all thresholds for filtering:
  neuroendocrine_enriched_clusters <- proportion_cells %>%
    left_join(cluster_avg_expr, by = c("seurat_clusters" = "cluster")) %>%
    filter(
      percentage > 5 &  # percentage of cells expressing > 5%
        log2fc > 0.25 # lfc > 0.25
    )
  
  # store results for  gene in the list
  nen_res_ARC[[gene]] <- neuroendocrine_enriched_clusters
  
  # average expression plot
  ggplot(avg_expression_melt, aes(x = as.factor(cluster), y = expression, fill = as.factor(cluster))) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic() +
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
          legend.position = "none") +
    labs(title = paste("Average expression of", gene),
         x = "Cluster", y = "Average expression") +
    scale_x_discrete(labels = 0:14) +
    scale_y_continuous(expand = c(0,0)) 
  
  ggsave(paste0(gene, "_avg_exp_ARC.png"), height = 4, width = 5)
  
  # percemtage of cells expressing plot
  ggplot(proportion_cells, aes(x = as.factor(seurat_clusters), y = percentage, fill = as.factor(seurat_clusters))) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic() +
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
          legend.position = "none") +
    labs(title = paste("Percentage of Cells Expressing", gene),
         x = "Cluster", y = "Percentage of Cells Expressing") +
    scale_x_discrete(labels = 0:14) +
    scale_y_continuous(expand = c(0, 0)) 
  
  ggsave(paste0(gene, "_per_exp_ARC.png"), height = 4, width = 5)
}

# --------------------------------------------------------------------------------
# unknown neuron cell type cluster
# --------------------------------------------------------------------------------

# subsetting neuronal clusters for subclustering
unknown_dev <- subset(nuclei_dev, idents = c("Unknown"))

# normalise and scale data as usual
unknown_dev <- NormalizeData(unknown_dev)
unknown_dev <- FindVariableFeatures(unknown_dev, 
                                    selection.method = "vst", 
                                    nfeatures = 2000)

top10_unknown_dev <- head(VariableFeatures(unknown_dev), 10)
plot_top10_unknown_dev <- VariableFeaturePlot(unknown_dev)

# remove peptides from highly variable genes
peptides <- c("CRH", "OXT", "AVP", "TRH", "SST", "GHRH", "POMC", "NPY", "AGRP", "KISS1", "TAC3", "GNRH1")
hvgs_unknown <- setdiff(VariableFeatures(unknown_dev), peptides)

# save variable genes
# ggsave("unknown_dev_top10_variable.png", height = 5, width = 7)

# continue with standard seurat workflow
unknown_dev <- ScaleData(unknown_dev,
                         features = hvgs_unknown,
                         vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
unknown_dev <- RunPCA(unknown_dev, features = hvgs_unknown)
ElbowPlot(unknown_dev)

# UMAP dimensionality reduction clusters
unknown_dev <- FindNeighbors(unknown_dev, dims = 1:20)
unknown_dev <- FindClusters(unknown_dev, resolution = 0.5)
unknown_dev <- RunUMAP(unknown_dev, dims = 1:20)

# visualising clusters
DimPlot(unknown_dev, reduction = "umap", group.by = "seurat_clusters", label = T)

# looking at cluster sizes
table(unknown_dev@meta.data$seurat_clusters) 
# ggsave("UMAP_unknown_unannotated.png", height = 5, width = 6)

# find differentially expressed genes from each cell type
unknown_dev_markers <- FindAllMarkers(unknown_dev, 
                                      only.unknowns = T, 
                                      min.pct = 0.25,
                                      min.diff.pct = 0.25,
                                      logfc.threshold = 0.25)

# get top 10 markers for each cluster
top10_markers <- function(x){
  unknown_dev_markers[unknown_dev_markers$cluster == x, ] %>% head(n = 20)}

# create df of top 10 markers for each cluster
top10_markers_unknown_df <- map_dfr(0:9, top10_markers)

# save as .csv file
setwd("/Volumes/ritd-ag-project-rd01lx-nhhan20/Michelle/analysed_data/neuroendocrine/csv/")
# write.csv(top10_markers_unknown_df, "unknown_dev_top10_markers.csv")

# change working directory back to plots
setwd("/Volumes/ritd-ag-project-rd01lx-nhhan20/Michelle/analysed_data/neuroendocrine/plots")

# identifying neuropeptidergic unknownpulations in the unknown
FeaturePlot(unknown_dev, features = c("CRH", "AVP", "OXT", "TRH", "SST", "PENK", "GNRH1"))
# ggsave("UMAP_unknown_neuropeptides.png", height = 25, width = 15, dpi = 300)

VlnPlot(unknown_dev, features = c("CRH", "AVP", "OXT", "TRH", "SST", "PENK", "GNRH1"))
# ggsave("Vln_unknown_neuropeptides.png", height = 15, width = 30)

# neuroendocrine co-markers
FeaturePlot(unknown_dev, features = c("NR3C1", "SCGN", "AGTR1", "TAC1", "TH", "ADARB2", "CALB1", "KCNMB4", "FOXP1",
                                      "THRB", "NFIB", "CALB2", "SFRP2", "OMP", "UCN3", "GHRH"))
# ggsave("UMAP_unknown_neuroendocrine.png", height = 20, width = 25, dpi = 300)

VlnPlot(unknown_dev, features = c("NR3C1", "SCGN", "PCSK1", "PCSK2", "AGTR1", "TAC1", "TH", "ADARB2", "CALB1", "KCNMB4", "FOXP1",
                                  "THRB", "NFIB", "CALB2", "SFRP2","OMP", "UCN3", "GHRH"))
# ggsave("vln_unknown_neuroendocrine.png", height = 20, width = 25)

# initialise list to store results in 
nen_res_unknown <- list()

# loop through each gene
for (gene in genes) {
  
  # pseudobulk expression data for gene per cluster
  avg_expression <- AverageExpression(unknown_dev,
                                      features = gene,
                                      return.seurat = F,
                                      assays = "RNA",
                                      group.by = "seurat_clusters")$RNA
  
  # convert to df and transpose to have clusters as rows
  avg_expression_df <- as.data.frame(t(avg_expression))
  
  # assign gene name as column name
  colnames(avg_expression_df) <- c(gene)
  
  # add cluster column
  avg_expression_df$cluster <- as.numeric(gsub("g", "", rownames(avg_expression_df)))
  
  # melt df for ggplot
  avg_expression_melt <- melt(avg_expression_df, id.vars = "cluster")
  
  # add column names for df
  colnames(avg_expression_melt) <- c("cluster", "gene", "expression")
  
  # calculate the proportion of cells expressing the gene (expression > 0)
  unknown_dev_data <- FetchData(unknown_dev, vars = c(gene, "seurat_clusters"))
  
  proportion_cells <- unknown_dev_data %>%
    group_by(seurat_clusters) %>%
    summarise(percentage = sum(!!sym(gene) > 0) / n() * 100)
  
  # convert "seurat_clusters" to numeric to match with "cluster" in avg_expression_melt
  proportion_cells$seurat_clusters <- as.numeric(as.character(proportion_cells$seurat_clusters))
  
  # calculate mean expression across all clusters
  mean_expr <- mean(avg_expression_melt$expression) 
  
  # calculate log2 fold change for each cluster relative to the global mean expression
  cluster_avg_expr <- avg_expression_melt %>%
    group_by(cluster) %>%
    summarise(avg_expr = mean(expression))
  
  # calculate log2 fold change (LFC) relative to the global mean expression
  cluster_avg_expr$log2fc <- log2(cluster_avg_expr$avg_expr / mean_expr)
  
  # apply all thresholds for filtering:
  neuroendocrine_enriched_clusters <- proportion_cells %>%
    left_join(cluster_avg_expr, by = c("seurat_clusters" = "cluster")) %>%
    filter(
      percentage > 5 &  # percentage of cells expressing > 5%
        log2fc > 0.25 # lfc > 0.25
    )
  
  # store results for  gene in the list
  nen_res_unknown[[gene]] <- neuroendocrine_enriched_clusters
  
  # average expression plot
  ggplot(avg_expression_melt, aes(x = as.factor(cluster), y = expression, fill = as.factor(cluster))) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic() +
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
          legend.position = "none") +
    labs(title = paste("Average expression of", gene),
         x = "Cluster", y = "Average expression") +
    scale_x_discrete(labels = 0:14) +
    scale_y_continuous(expand = c(0,0)) 
  
  ggsave(paste0(gene, "_avg_exp_unknown.png"), height = 4, width = 5)
  
  # percemtage of cells expressing plot
  ggplot(proportion_cells, aes(x = as.factor(seurat_clusters), y = percentage, fill = as.factor(seurat_clusters))) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic() +
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
          legend.position = "none") +
    labs(title = paste("Percentage of Cells Expressing", gene),
         x = "Cluster", y = "Percentage of Cells Expressing") +
    scale_x_discrete(labels = 0:14) +
    scale_y_continuous(expand = c(0, 0)) 
  
  ggsave(paste0(gene, "_per_exp_.png"), height = 4, width = 5)
}

# --------------------------------------------------------------------------------
# labelling neuroendocrine populations
# --------------------------------------------------------------------------------

# create label for NEN
nuclei_dev$NEN_type <- "non_NEN" # label every cell as non_NEN before rewriting NEN labels

# add labels for PVN NENs
nuclei_dev$NEN_type[colnames(nuclei_dev) %in% 
                              colnames(PVN_dev)[PVN_dev$seurat_clusters %in% c(5, 7, 11, 13, 14)]] <- "PVN_NEN"

# add labels for other NENs
nuclei_dev$NEN_type[colnames(nuclei_dev) %in% 
                              colnames(PO_dev)[PO_dev$seurat_clusters %in% c(5, 11)]] <- "other_NEN"

nuclei_dev$NEN_type[colnames(nuclei_dev) %in% 
                              colnames(ARC_dev)[ARC_dev$seurat_clusters %in% 5]] <- "other_NEN"

nuclei_dev$NEN_type[colnames(nuclei_dev) %in% 
                              colnames(AH_dev)[AH_dev$seurat_clusters %in% c(5, 11, 12)]] <- "other_NEN"

# add labels for NEN subtypes
nuclei_dev$NEN_subtype <- "non_NEN" 

nuclei_dev$NEN_subtype[colnames(nuclei_dev) %in% 
                                 colnames(PVN_dev)[PVN_dev$seurat_clusters %in% 5]] <- "TRH_1"

nuclei_dev$NEN_subtype[colnames(nuclei_dev) %in% 
                                 colnames(PVN_dev)[PVN_dev$seurat_clusters %in% 7]] <- "CRH/AVP"

nuclei_dev$NEN_subtype[colnames(nuclei_dev) %in% 
                                 colnames(PVN_dev)[PVN_dev$seurat_clusters %in% 11]] <- "AVP"

nuclei_dev$NEN_subtype[colnames(nuclei_dev) %in% 
                                 colnames(PVN_dev)[PVN_dev$seurat_clusters %in% 13]] <- "AVP/OXT"

nuclei_dev$NEN_subtype[colnames(nuclei_dev) %in% 
                                 colnames(PVN_dev)[PVN_dev$seurat_clusters %in% 14]] <- "TRH_2"

nuclei_dev$NEN_subtype[colnames(nuclei_dev) %in% 
                                 colnames(PO_dev)[PO_dev$seurat_clusters %in% c(5, 11)]] <- "GNRH1"

nuclei_dev$NEN_subtype[colnames(nuclei_dev) %in% 
                                 colnames(ARC_dev)[ARC_dev$seurat_clusters %in% 5]] <- "GHRH"

nuclei_dev$NEN_subtype[colnames(nuclei_dev) %in% 
                                 colnames(AH_dev)[AH_dev$seurat_clusters %in% 5]] <- "TRH_3"

nuclei_dev$NEN_subtype[colnames(nuclei_dev) %in% 
                                 colnames(AH_dev)[AH_dev$seurat_clusters %in% c(11, 12)]] <- "TRH_4"

# create stacked bar plot per region of NENs vs non-NENs
NEN_region_summary <- nuclei_dev@meta.data %>%
  filter(cell_type != "Unknown") %>%
  group_by(cell_type, NEN_type) %>%
  summarise(cell_count = n()) %>%
  ungroup()

NEN_region_summary$NEN_type <- ifelse(NEN_region_summary$NEN_type == "non_NEN", "non_NEN", "NEN") # keep NEN and non-NEN labels only for stacked bar chart

NEN_region_summary$cell_type <- factor(NEN_region_summary$cell_type, levels = c("Prethalamus", "Thalamus", "SMN", "SCN", "LH",
                                                                                "DMH", "VMH", "PH", "MN", "TMN", "TH",
                                                                                "ARC","PO", "AH", "PVN"))

ggplot(NEN_region_summary, aes(x = cell_type, y = cell_count, fill = NEN_type)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(title = "Proportion of NENs vs non-NENs per region",
       X = "Region", y = "Cell count") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0)) +
  coord_flip()

# ggsave("NEN_per_region.png", height = 4, width = 5, dpi = 300)

# create a bar chart for global proportions of NENs vs non-NENs
NEN_summary <- nuclei_dev@meta.data %>%
  filter(cell_type != "Unknown") %>%
  group_by(NEN_type) %>%
  summarise(cell_count = n())

ggplot(NEN_summary, aes(x = NEN_type, y = cell_count, fill = NEN_type)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(title = "Proportion of NENs vs non-NENs per region",
       X = "Region", y = "Cell count") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0)) +
  coord_flip()

# ggsave("NEN_proportion_global.png", height = 4, width = 5, dpi = 300)

# create pie chart for global prportions of NENs and non-NENs
NEN_summary_pie <- NEN_summary %>%
  mutate(percentage = cell_count / sum(cell_count) * 100)

ggplot(NEN_summary_pie, aes(x = "", y = percentage, fill = NEN_type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void()

# ggsave("NEN_pie_proportions.png", height = 5, width = 5, dpi = 300)

# plot violin plot of comarkers for PVN clusters
PVN_NEN <- subset(PVN_dev, idents = c(5, 7, 11, 13, 14))

PVN_NEN <- RenameIdents(PVN_NEN,
                        c("5" = "TRH_1", 
                          "7" = "CRH/AVP",
                          "11" = "AVP", 
                          "13" = "AVP/OXT",
                          "14" = "TRH_2"))

comarkers <- c("NR3C1", "TH", "KCNMB4", "FOXP1", "TAC1", "NFIB", "CALB1")

for (gene in comarkers) {
  # plot violin plot amongst PVN NEN clusters for each gene
  VlnPlot(PVN_NEN, features = gene)
  
  # save plot
  ggsave(paste0("PVN_NEN_", gene, ".png"), height = 5, width = 6, dpi = 300)
  
}

# plot violin plot of comarkers for PO clusters
PO_NEN <- subset(PO_dev, idents = c(5, 11))

PO_NEN <- RenameIdents(PO_NEN,
                       c("5" = "GNRH1", 
                         "11" = "GNRH1"))

comarkers <- c("SP4")

for (gene in comarkers) {
  # plot violin plot amongst PO NEN clusters for each gene
  VlnPlot(PO_NEN, features = gene)
  
  # save plot
  ggsave(paste0("PO_NEN_", gene, ".png"), height = 5, width = 6, dpi = 300)
  
}

# plot violin plot of comarkers for ARC clusters
ARC_NEN <- subset(ARC_dev, idents = c(5))

ARC_NEN <- RenameIdents(ARC_NEN,
                        "5" = "GNRH")

comarkers <- c("FOXP2", "ADARB2")

for (gene in comarkers) {
  # plot violin plot amongst ARC NEN clusters for each gene
  VlnPlot(ARC_NEN, features = gene)
  
  # save plot
  ggsave(paste0("ARC_NEN_", gene, ".png"), height = 5, width = 6, dpi = 300)
  
}

# plot violin plot of comarkers for AH clusters
AH_NEN <- subset(AH_dev, idents = c(5, 11, 12))

AH_NEN <- RenameIdents(AH_NEN,
                       c("5" = "TRH_1", 
                         "11" = "TRH_2", 
                         "12" = "TRH_2"))

comarkers <- c("TAC1", "GHRH")

for (gene in comarkers) {
  # plot violin plot amongst AH NEN clusters for each gene
  VlnPlot(AH_NEN, features = gene)
  
  # save plot
  ggsave(paste0("AH_NEN_", gene, ".png"), height = 5, width = 6, dpi = 300)
  
}

# create dot plot of co-marker expression for NEN identification
nuclei_dev <- SetIdent(nuclei_dev, value = "NEN_type")

NEN_clusters <- subset(nuclei_dev, idents = c("PVN_NEN", "other_NEN")) # subset NEN clusters
NEN_clusters <- SetIdent(NEN_clusters, value = "NEN_subtype")

NEN_clusters$NEN_subtype <- factor(
  NEN_clusters$NEN_subtype, 
  levels = c("CRH/AVP", "AVP", "AVP/OXT", "TRH_1", "TRH_2", "TRH_3", "TRH_4", "GHRH", "GNRH1")  # reorder subtype populations
)

markers <- c("SP4", "GNRH1", "ADARB2", "FOXP2", "GHRH", "TAC1", "TRH", "FOXP1", "OXT", "TH", "KCNMB4", "CALB1", "AVP",
             "NR3C1", "CRH")

DotPlot(NEN_clusters, features = markers, group.by = "NEN_subtype") +
  scale_size(range = c(0, 8)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Marker expression across NEN clusters") + 
  coord_flip()

# ggsave("NEN_marker_dotplot.png", height = 5, width = 6, dpi = 300)

# create a heatmap for marker identification 
expr_matrix <- as.matrix(NEN_clusters[["RNA"]]@data[markers, ]) # extract expression matrix (scaled counts)
expr_matrix_scaled <- t(scale(t(expr_matrix))) # scale per gene

cluster_annotation <- data.frame(
  Cluster = Idents(NEN_clusters)
)
rownames(cluster_annotation) <- colnames(NEN_clusters)

ha <- HeatmapAnnotation(df = cluster_annotation, 
                        col = list(Cluster = c("CRH/AVP" = "#E41A1C",
                                               "AVP" = "#377EB8",
                                               "AVP/OXT" = "#4DAF4A",
                                               "TRH_1" = "#984EA3",
                                               "TRH_2" = "#FF7F00",
                                               "TRH_3" = "#FFFF33",
                                               "TRH_4" = "#a65628",
                                               "GHRH" = "#f781bf",
                                               "GNRH1" = "#999999")))
# Draw heatmap
png("NEN_heatmap.png", height = 1500, width = 2000, res = 300)
Heatmap(expr_matrix_scaled,
        name = "Scaled Expression",
        top_annotation = ha,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = TRUE,
        row_names_side = "left",
        col = colorRamp2(c(-1, 0, 2), c("white","lightblue","blue")),
        column_split = NEN_clusters$NEN_subtype,
        column_title = NULL)
dev.off()
