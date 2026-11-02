# load required packages
library(here)
library(tidyverse)
library(Seurat)
library(harmony)

# set working directory
setwd(here("plots"))

# run UMAP
seurat_dev <- FindNeighbors(seurat_dev, reduction = "harmony", dims = 1:40)
seurat_dev <- FindClusters(seurat_dev, resolution = 0.4)
seurat_dev <- RunUMAP(seurat_dev, reduction = "harmony", dims = 1:40)

# visualise UUMAP
DimPlot(seurat_dev, group.by = "seurat_clusters", reduction = "umap", label = TRUE) 
ggsave("umap_major_unannotated.png", height = 7, width = 8)

# cluster sizes
table(seurat_dev@meta.data$seurat_clusters) 

# neuroepithelial markers
FeaturePlot(seurat_dev, features = c("SOX2", "SOX1", "SOX10", "PAX6", "VIM", "TJP1", "OCLN", "NES", "HES1", "NOTCH1", "CDH1",
                                     "HMGA2", "ARHGAP28"))
ggsave("umap_neuroepithelial_markers.png", height = 25, width = 30, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("SOX2", "SOX1", "SOX10", "PAX6", "VIM", "TJP1", "OCLN", "NES", "HES1", "NOTCH1", "CDH1",
                     "HMGA2", "ARHGAP28"), 
        pt.size = 0)
ggsave("violin_neuroepithelial_markers.png", height = 15, width = 35, dpi = 300)

# radial glia markers
FeaturePlot(seurat_dev, 
            features = c("SOX2", "PAX6", "VIM", "NES", "HES1", "HES5", "CDH2", "TNC", "GFAP", "FABP7", "SLC1A3", "HMGA2",
                         "COL2A1", "TTYH1", "MKI67", "TOP2A"))
ggsave("umap_radialglia_markers.png", height = 25, width = 30, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("SOX2", "PAX6", "VIM", "NES", "HES1", "HES5", "CDH2", "TNC", "GFAP", "FABP7", "SLC1A3", "HMGA2",
                     "COL2A1", "TTYH1", "MKI67", "TOP2A"), 
        pt.size = 0)
ggsave("violin_radialglia_markers.png", height = 15, width = 35, dpi = 300)

# neural progenitor cell (NPC) / intermediate progenitor markers
FeaturePlot(seurat_dev, 
            features = c("ASCL1", "NEUROG2", "NEUROD1", "EOMES", "DCX", "FABP7", "NES", "NBPF10", "NBPF15", "STMN2", "PLAGL1", 
                         "SFTA3", "EGFR", "SYT1", "NHLH1", "NHLH2"))
ggsave("umap_npc_markers.png", height = 25, width = 25, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("ASCL1", "NEUROG2", "NEUROD1", "EOMES", "DCX", "FABP7", "NES", "NBPF10", "NBPF15", "STMN2", "PLAGL1", 
                     "SFTA3", "EGFR", "SYT1", "NHLH1", "NHLH2"), 
        pt.size = 0)
ggsave("violin_npc_markers.png", height = 15, width = 35, dpi = 300)

# neuron markers
FeaturePlot(seurat_dev, 
            features = c("SYT1", "DCX", "STMN2", "TUBB3", "MAP2", "RBFOX3", "SYP", "GAP43", "NEFL", "NEFM", "SNAP25"))
ggsave("umap_neuron_markers.png", height = 18, width = 25, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("SYT1", "DCX", "STMN2", "TUBB3", "MAP2", "RBFOX3", "SYP", "GAP43", "NEFL", "NEFM", "SNAP25"), 
        pt.size = 0)
ggsave("violin_neuron_markers.png", height = 10, width = 30, dpi = 300)

# astrocyte markers
FeaturePlot(seurat_dev, 
            features = c("GFAP", "SLC1A2", "SLC1A3", "GLUL", "S100B", "ALDH1L1", "AQP4", "GJA1", "SOX9", "FAM107A", "HOPX", 
                         "AGT", "ALDOC"))
ggsave("umap_astrocyte_markers.png", height = 25, width = 25, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("GFAP", "SLC1A2", "SLC1A3", "GLUL", "S100B", "ALDH1L1", "AQP4", "GJA1", "SOX9", "FAM107A", "HOPX", 
                     "AGT", "ALDOC"), 
        pt.size = 0)
ggsave("violin_astrocyte_markers.png", height = 15, width = 30, dpi = 300)

# oligodendrocyte precursor cell (OPC) markers
FeaturePlot(seurat_dev, 
            features = c("PDGFRA", "CSPG4", "OLIG1", "OLIG2", "SOX10", "NKX2-2", "EGFR", "GPR75-ASB3", "APOD"))
ggsave("umap_OPC_markers.png", height = 10, width = 12, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("PDGFRA", "CSPG4", "OLIG1", "OLIG2", "SOX10", "NKX2-2", "EGFR", "GPR75-ASB3", "APOD"), 
        pt.size = 0)
ggsave("violin_OPC_markers.png", height = 10, width = 25, dpi = 300)

# oligodendrocyte markers
FeaturePlot(seurat_dev, 
            features = c("MBP", "MOG", "MOBP", "MAG", "PLP1", "CNP", "CLDN11", "OLIG1", "OLIG2", "BCAS1"))
ggsave("umap_oligo_markers.png", height = 18, width = 25, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("MBP", "MOG", "MOBP", "MAG", "PLP1", "CNP", "CLDN11", "OLIG1", "OLIG2", "BCAS1"), 
        pt.size = 0)
ggsave("violin_oligo_markers.png", height = 10, width = 30, dpi = 300)

# microglia
FeaturePlot(seurat_dev, 
            features = c("C1QA", "C1QB", "C1QC", "TYROBP", "AIF1", "P2RY12", "TMEM119", "CX3CR1", "CSF1R", "ITGAM"))
ggsave("umap_microglia_markers.png", height = 18, width = 25, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("C1QA", "C1QB", "C1QC", "TYROBP", "AIF1", "P2RY12", "TMEM119", "CX3CR1", "CSF1R", "ITGAM"), 
        pt.size = 0)
ggsave("violin_microglia_markers.png", height = 10, width = 30, dpi = 300)

# tanycyte markers
FeaturePlot(seurat_dev, 
            features = c("RAX", "DIO2", "FGF10", "VIM", "SLC1A3", "CRYM", "SLC16A2"))
ggsave("umap_tanycyte_markers.png", height = 10, width = 12, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("RAX", "DIO2", "FGF10", "VIM", "SLC1A3", "CRYM", "SLC16A2"), 
        pt.size = 0)
ggsave("violin_tanycyte_markers.png", height = 10, width = 25, dpi = 300)

# ependymal cell markers
FeaturePlot(seurat_dev, 
            features = c("FOXJ1", "DNAI1", "DNAH5", "DNAH9", "CCDC74A", "CCDC74B", "LRTOMT", "CCDC153"))
ggsave("umap_ependymal_markers.png", height = 10, width = 12, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("FOXJ1", "DNAI1", "DNAH5", "DNAH9", "CCDC74A", "CCDC74B", "LRTOMT", "CCDC153"), 
        pt.size = 0)
ggsave("violin_ependymal_markers.png", height = 10, width = 30, dpi = 300)

# endothelial markers
FeaturePlot(seurat_dev, 
            features = c("CLDN5", "PECAM1", "VWF", "FLT1", "CD34", "KDR", "EMCN", "SLC38A5", "ITM2A"))
ggsave("umap_endothelial_markers.png", height = 18, width = 25, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("CLDN5", "PECAM1", "VWF", "FLT1", "CD34", "KDR", "EMCN", "SLC38A5", "ITM2A"), 
        pt.size = 0)
ggsave("violin_endothelial_markers.png", height = 10, width = 25, dpi = 300)

# vascular leptomenegial cell (VLMC) markers
FeaturePlot(seurat_dev, 
            features = c("PTGDS", "COL1A1", "COL1A2", "COL3A1", "LUM", "DCN"))
ggsave("umap_vlmc_markers.png", height = 10, width = 10, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("PTGDS", "COL1A1", "COL1A2", "COL3A1", "LUM", "DCN"), 
        pt.size = 0)
ggsave("violin_vlmc_markers.png", height = 8, width = 20, dpi = 300)

# mural markers
FeaturePlot(seurat_dev, 
            features = c("PDGFRB", "MYH11", "CNN1", "MCAM", "DES", "ACTA2", "TAGLN", "KCNJ8", "KCNJ6", "ABCC9", "RGS5"))
ggsave("umap_mural_markers.png", height = 18, width = 25, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("PDGFRB", "MYH11", "CNN1", "MCAM", "DES", "ACTA2", "TAGLN", "KCNJ8", "KCNJ6", "ABCC9", "RGS5"), 
        pt.size = 0)
ggsave("violin_mural_markers.png", height = 10, width = 30, dpi = 300)

# find gene markers for each cluster
major_markers <- FindAllMarkers(seurat_dev, 
                                only.pos = T, 
                                min.pct = 0.25,
                                min.diff.pct = 0.25,
                                logfc.threshold = 0.25)

# get top 10 markers for each cluster
top10_markers <- function(x){
  major_markers[major_markers$cluster == x, ] %>% head(n = 10)
  }

top10_markers_df <- map_dfr(0:41, top10_markers)

# save as .csv file
write.csv(top10_markers_df, here("csv", "major_top10_markers.csv"))

# renaming cluster to cell types
new_cluster_ids <- c("Neuron", "Astrocyte", "OPC", "Neuron", "Neuron", "Astrocyte", "RG", "Microglia", "Neuron", "RG",
                     "RG", "IPC", "Neuron", "Neuroepithelial", "Oligodendrocyte", "Pericyte", "Neuron", "Ependymal", "Endothelial", 
                     "Cycling progenitor", "Ependymal", "Cycling progenitor", "Neuroepithelial", "Ependymal", "VLMC", "Cycling progenitor", 
                     "OPC", "Microglia", "Neuron", "Astrocyte", "Astrocyte", "Neuron", "Cycling progenitor", "Neuron", "Macrophage", "Neuron", 
                     "Neuron", "Neuron", "Neuron", "Cycling progenitor", "Cycling progenitor", "Tanycyte")

names(new_cluster_ids) <- levels(seurat_dev)
seurat_dev <- RenameIdents(seurat_dev, new_cluster_ids)
seurat_dev@meta.data$"cell_type" <- as.factor(seurat_dev@active.ident)
seurat_dev$cell_type <- factor(seurat_dev$cell_type, levels = sort(unique(as.character(seurat_dev$cell_type))))

# visualise annotated UMAP
seurat_dev$cell_type <- factor(seurat_dev$cell_type,
                               levels = c("Neuroepithelial", "RG", "IPC", "Cycling progenitor", "Neuron",
                                          "Astrocyte", "OPC", "Oligodendrocyte", "Microglia", "Macrophage", "Ependymal", "Tanycyte",
                                          "Endothelial", "Pericyte", "VLMC"))

DimPlot(seurat_dev, reduction = "umap", group.by = "cell_type", label = FALSE, pt.size = 1.5) +
  scale_colour_manual(
    values = c("Astrocyte" = "#D78FA9", "Cycling progenitor" = "#8E3F54", "Endothelial" = "#755333", "Ependymal" = "#829B73",
               "IPC" = "#D5A928", "Macrophage" = "#B85F78", "Microglia" = "#765996", "Neuroepithelial" = "#DD805E", 
               "Neuron" = "#A995C9", "Oligodendrocyte" = "#B2A52E", "OPC" = "#5FA9A2", "Pericyte" = "#3F8C5A", "RG" = "#78A8D1",
               "Tanycyte" = "#587F9E", "VLMC" = "#B85F5A"),
    labels = function(x) sub(" progenitor$", "\nprogenitor", x)) +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  theme_void() +
  theme(plot.title = element_blank(), 
        legend.text = element_text(size = 7),
        legend.key.height = unit(4, "mm"),
        legend.box.spacing = unit(0, "pt"))

ggsave("umap_major_annotated.tiff", height = 2.5, width = 3.2, dpi = 600, compression = "lzw")

# plot dotplot of markers
DotPlot(seurat_dev, assay = "RNA", group.by = "cell_type",
        features = c("LUM", "COL1A1", "PTGDS",             # VLMC
                     "ABCC9", "KCNJ8", "PDGFRB",            # pericyte
                     "PECAM1", "FLT1", "CLDN5",            # endothelial
                     "CRYM",                               # tanycyte
                     "FOXJ1", "DNAI1", "DNAH9",            # ependymal
                     "CD163", "LYVE1", "MRC1",             # macrophage
                     "TMEM119", "CX3CR1", "P2RY12",        # microglia
                     "MOG", "MBP", "BCAS1",                # oligodendrocyte
                     "PCDH15", "PDGFRA", "OLIG1",          # OPC
                     "GFAP", "SLC1A2", "AQP4",             # astrocyte
                     "GAP43", "TUBB3", "STMN2",            # neuron
                     "PCNA", "TOP2A", "MKI67",             # cycling progenitor
                     "HES6", "EGFR", "ASCL1",              # intermediate progenitor
                     "COL2A1", "HES1", "VIM",              # radial glia
                     "PAX6", "HMGA2", "CDH2"               # neuroepithelial
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
        legend.key.height = unit(5, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.box.spacing = unit(0, "pt"),
        panel.grid = element_blank(),
        panel.heights = unit(43 * 0.3, "cm"),
        panel.widths = unit(length(unique(seurat_dev$cell_type)) * 0.4, "cm"))

ggsave("major_marker_dotplot.tiff", height = 6, width = 3.5, dpi = 600, compression = "lzw")

# remove objects
rm(new_cluster_id)