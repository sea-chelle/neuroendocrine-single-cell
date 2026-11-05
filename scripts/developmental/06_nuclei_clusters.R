# Load required packages
library(Seurat)
library(tidyverse)
library(here)

# set working directory to figures
setwd(here("plots"))

# exclude immature neurons for hypothalamic nuclei identification
nuclei_dev <- subset(neuron_dev, subset = cell_type != "Immature")

# normalise data and find variable features
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
nuclei_dev <- FindClusters(nuclei_dev, resolution = 0.5)
nuclei_dev <- RunUMAP(nuclei_dev, dims = 1:40)

# visualising clusters
DimPlot(nuclei_dev, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
ggsave("UMAP_neuroendocrine_unannotated.png", height = 7, width = 8)

# looking at cluster sizes
table(nuclei_dev@meta.data$seurat_clusters) 

# AH markers
FeaturePlot(nuclei_dev, features = c("SIX3", "DLX1", "DLX2", "DLX5", "DLX6", "ACHE"))
ggsave("umap_AH_marker.png", height = 20, width = 15, dpi = 300)

VlnPlot(nuclei_dev, features = c("SIX3", "DLX1", "DLX2", "DLX5", "DLX6", "ACHE"), pt.size = 0)
ggsave("violin_AH_marker.png", height = 10, width = 30, dpi = 300)

# PVN markers
FeaturePlot(nuclei_dev, features = c("SIM1", "OTP", "POU3F2"))
ggsave("umap_PVN_marker.png", height = 10, width = 20, dpi = 300)

VlnPlot(nuclei_dev, features = c("SIM1", "OTP", "POU3F2"), pt.size = 0)
ggsave("violin_PVN_marker.png", height = 5, width = 30, dpi = 300)

# POA markers
FeaturePlot(nuclei_dev, features = c("ARX", "NFIX", "GAL"))
ggsave("umap_POA_marker.png", height = 10, width = 20, dpi = 300)

VlnPlot(nuclei_dev, features = c("ARX", "NFIX", "GAL"), pt.size = 0)
ggsave("violin_POA_marker.png", height = 5, width = 30, dpi = 300)

# SCN markers
FeaturePlot(nuclei_dev, features = c("RGS16", "RELN", "ONECUT1", "ONECUT2", "SIX6", "LHX1", "DLX2", "SIX3"))
ggsave("umap_SCN_marker.png", height = 20, width = 30, dpi = 300)

VlnPlot(nuclei_dev, features = c("RGS16", "RELN", "ONECUT1", "ONECUT2", "SIX6", "LHX1", "DLX2", "SIX3"), pt.size = 0)
ggsave("violin_SCN_marker.png", height = 15, width = 35, dpi = 300)

# TH markers
FeaturePlot(nuclei_dev, features = c("ISL1", "MEIS2", "SST", "DLK1", "TBX3"))
ggsave("umap_TH_marker.png", height = 20, width = 15, dpi = 300)

VlnPlot(nuclei_dev, features = c("ISL1", "MEIS2", "SST", "DLK1", "TBX3"), pt.size = 0)
ggsave("violin_TH_marker.png", height = 10, width = 30, dpi = 300)

# ARC markers
FeaturePlot(nuclei_dev, features = c("TBX3", "ISL1", "AGRP", "POMC", "GHRH", "KISS1", "TAC3", "NPY"))
ggsave("umap_ARC_marker.png", height = 20, width = 25, dpi = 300)

VlnPlot(nuclei_dev, features = c("TBX3", "ISL1", "AGRP", "POMC", "GHRH", "KISS1", "TAC3", "NPY"), pt.size = 0)
ggsave("violin_ARC_marker.png", height = 15, width = 30, dpi = 300)

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

# LH markers
FeaturePlot(nuclei_dev, features = c("LHX9", "PDYN", "HCRT"))
ggsave("umap_LH_marker.png", height = 20, width = 20, dpi = 300)

VlnPlot(nuclei_dev, features = c("LHX9", "PDYN", "HCRT"), pt.size = 0)
ggsave("violin_LH_marker.png", height = 5, width = 30, dpi = 300)

# TMN markers
FeaturePlot(nuclei_dev, features = c("HDC", "SLC18A2", "CBLN4", "TBX3"))
ggsave("umap_TMN_marker.png", height = 20, width = 20, dpi = 300)

VlnPlot(nuclei_dev, features = c("HDC", "SLC18A2", "CBLN4", "TBX3"), pt.size = 0)
ggsave("violin_TMN_marker.png", height = 10, width = 30, dpi = 300)

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

# MN markers
FeaturePlot(nuclei_dev, features = c("FOXB1", "RPRM", "CALB1", "LHX5", "PITX2"))
ggsave("umap_MN_marker.png", height = 20, width = 15, dpi = 300)

VlnPlot(nuclei_dev, features = c("FOXB1", "RPRM", "CALB1", "LHX5", "PITX2"), pt.size = 0)
ggsave("violin_MN_marker.png", height = 10, width = 30, dpi = 300)

# SMN markers
FeaturePlot(nuclei_dev, features = c("LMX1A", "PITX2", "CALB2", "BARHL1", "BARHL2", "NEUROD2"))
ggsave("umap_SMN_marker.png", height = 20, width = 15, dpi = 300)

VlnPlot(nuclei_dev, features = c("LMX1A", "PITX2", "CALB2", "BARHL1", "BARHL2", "NEUROD2"), pt.size = 0)
ggsave("violin_SMN_marker.png", height = 10, width = 30, dpi = 300)

# TE markers
FeaturePlot(nuclei_dev, features = c("TBR1", "EOMES", "MEIS2"))
ggsave("umap_TE_marker.png", height = 20, width = 20, dpi = 300)

VlnPlot(nuclei_dev, features = c("TBR1", "EOMES", "MEIS2"), pt.size = 0)
ggsave("violin_TE_marker.png", height = 5, width = 30, dpi = 300)

# ZI markers
FeaturePlot(nuclei_dev, features = c("MEIS2", "GABRA5"))
ggsave("umap_ZI_marker.png", height = 10, width = 20, dpi = 300)

VlnPlot(nuclei_dev, features = c("MEIS2"), pt.size = 0)
ggsave("violin_ZI_marker.png", height = 5, width = 30, dpi = 300)

# Prethalamus markers
FeaturePlot(nuclei_dev, features = c("SP9", "ARX"))
ggsave("umap_prethalamus_marker.png", height = 10, width = 20, dpi = 300)

VlnPlot(nuclei_dev, features = c("SP9", "ARX"), pt.size = 0)
ggsave("violin_prethalamus_marker.png", height = 5, width = 30, dpi = 300)

# thalamus markers
FeaturePlot(nuclei_dev, features = c("GBX2", "TCF7L2", "GATA2", "GATA3"))
ggsave("umap_thalamus_marker.png", height = 20, width = 20, dpi = 300)

VlnPlot(nuclei_dev, features = c("GBX2", "TCF7L2", "GATA2", "GATA3"), pt.size = 0)
ggsave("violin_thalamus_marker.png", height = 10, width = 30, dpi = 300)

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
top10_markers_neuroendocrine_df <- map_dfr(0:34, top10_markers)

# save as .csv file
write.csv(top10_markers_neuroendocrine_df, here("csv", "nuclei_dev_top10_markers.csv"))

# renaming cluster numbers to cell types
new_cluster_ids_nuclei <- c("Thalamus", "AH", "Tuberal", "Unresolved", "Unresolved", "Unresolved", "TMN", "PVN", "POA", "VMH",
                            "PMN", "Unresolved", "Prethalamus", "ARC", "PH", "Thalamus", "Unresolved", "ARC", "Tuberal", "Unresolved",
                            "TMN", "PVN", "LH", "ARC", "Telencephalic", "SCN", "PVN", "Prethalamus", "PH", "Thalamus",
                            "SMN", "Thalamus", "Telencephalic", "MB", "Telencephalic")

names(new_cluster_ids_nuclei) <- levels(nuclei_dev)
nuclei_dev <- RenameIdents(nuclei_dev, new_cluster_ids_nuclei)
nuclei_dev@meta.data$"cell_type" <- as.factor(nuclei_dev@active.ident)
nuclei_dev$cell_type <- factor(nuclei_dev$cell_type, levels = sort(unique(as.character(nuclei_dev$cell_type))))

DimPlot(nuclei_dev, reduction = "umap", group.by = "cell_type", label = FALSE, pt.size = 0.01)

# visualise annotated UMAP
DimPlot(nuclei_dev, reduction = "umap", group.by = "cell_type", label = FALSE, pt.size = 0.01) +
  scale_colour_manual(values = c("AH" = "#D78FA9", "ARC" = "#5FA9A2", "LH" = "#8E3F54", "MB" = "#B85F78", "PH" = "#765996", 
                                 "PMN" = "#9BCB9A", "POA" = "#78A8D1", "Prethalamus" = "#587F9E",  "PVN" = "#A995C9", "SCN" = "#DD805E", 
                                 "SMN" = "#829B73", "Telencephalic" = "#B85F5A",  "Thalamus" = "#3F8C5A", "TMN" = "#D5A928", 
                                 "Tuberal" = "#B2A52E", "Unresolved" = "grey80", "VMH" = "#755333")) +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  theme_void() +
  theme(plot.title = element_blank(), 
        legend.text = element_text(size = 7),
        legend.key.height = unit(4, "mm"),
        legend.box.spacing = unit(0, "pt"))

ggsave("umap_nuclei_annotated.tiff", height = 2.5, width = 3.1, dpi = 600, compression = "lzw")

DimPlot(nuclei_dev, reduction = "umap", group.by = "cell_type", label = FALSE, pt.size = 0.0001) +
  scale_colour_manual(values = c("AH" = "#DD805E", "ARC" = "#D78FA9", "LH" = "#D5A928", "MB" = "#B85F78", "PH" = "#B85F5A", 
                                 "PMN" = "#3F8C5A", "POA" = "#78A8D1", "Prethalamus" = "#587F9E",  "PVN" = "#A995C9", "SCN" = "#B2A52E", 
                                 "SMN" = "#9BCB9A", "Telencephalic" = "#755333",  "Thalamus" = "#5FA9A2", "TMN" = "#765996", 
                                 "Tuberal" = "#829B73", "Unresolved" = "grey80", "VMH" = "#8E3F54")) +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  theme_void() +
  theme(plot.title = element_blank(), 
        legend.text = element_text(size = 7),
        legend.key.height = unit(4, "mm"),
        legend.box.spacing = unit(0, "pt"))
ggsave("umap_nuclei_annotated.tiff", height = 2.8, width = 3, dpi = 600, compression = "lzw")

table(nuclei_dev@meta.data$cell_type) # looking at cluster sizes

# set order
nuclei_dev$cell_type <- factor(nuclei_dev$cell_type,
                               levels = c("AH", "PVN", "POA", "SCN", "Tuberal", "ARC", "TMN", "VMH", "LH", "PH", 
                                          "MB", "PMN", "SMN", "Prethalamus", "Thalamus", "Telencephalic", "Unresolved"))

# plot dotplot of markers
DotPlot(nuclei_dev, assay = "RNA", group.by = "cell_type",
        features = c("BCL11B", "SOX6", "MEF2C", "FOXG1", # Telencephalic
                     "GBX2", "TCF7L2",  # Thalamus
                     "MEIS2", "SP9",  # Prethalamus
                     "BARHL2", "LMX1A",  # SMN
                     "PITX2", "ADCYAP1", # PMN
                     "RPRM", "FOXB1", # MB
                     "CRABP1", "IRX2",  # PH
                     "PDYN", "HCRT",  # LH
                     "NR5A1", "FEZF1",  # VMH
                     "SLC18A2", "HDC", # TMN
                     "KISS1", "POMC", "GHRH", # ARC
                     "TBX3", "ISL1", # Tuberal
                     "SIX3", "RELN", # SCN
                     "NFIX",  # POA
                     "OTP", "SIM1", # PVN
                     "DLX2", "ARX"  # AH
                     ),
        dot.scale = 3, cols = c("grey100", "#7656A8")) +
  coord_flip() +
  guides(colour = guide_colourbar(title = "Average\nexpression"),
         size = guide_legend(title = "Percent\nexpressed")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7, face = "italic"),
        axis.title = element_blank(),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(3, "mm"),
        legend.key.height = unit(4, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.box.spacing = unit(0, "pt"),
        panel.grid = element_blank(),
        panel.heights = unit(34 * 0.3, "cm"),
        panel.widths = unit(length(unique(seurat_dev$cell_type)) * 0.4, "cm"))

ggsave("nuclei_marker_dotplot.tiff", height = 4.6, width = 3.5, dpi = 600, compression = "lzw")