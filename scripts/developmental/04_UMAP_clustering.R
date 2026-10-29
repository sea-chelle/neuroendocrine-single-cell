# load required packages
library(here)
library(Seurat)
library(harmony)

# set working directory
setwd(here("plots"))

# run UMAP
seurat_dev <- FindNeighbors(seurat_dev, reduction = "harmony", dims = 1:50)
seurat_dev <- FindClusters(seurat_dev, resolution = 0.5)
seurat_dev <- RunUMAP(seurat_dev, reduction = "harmony", dims = 1:50)

# visualise UUMAP
DimPlot(seurat_dev, group.by = "seurat_clusters", reduction = "umap", label = TRUE) 
table(seurat_dev@meta.data$seurat_clusters) # look at cluster sizes

ggsave("umap_major_unannotated.png", height = 7, width = 8)

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

top10_markers_df <- map_dfr(0:34, top10_markers)

# save as .csv file
setwd(here("csv"))
write.csv(top10_markers_df, "major_top10_markers.csv")

# neuroepithelial
FeaturePlot(seurat_dev, 
            features = c("NES", "VIM", "PAX6", "SOX2", "SOX10", "NOTCH1", "HES1", "OCLN", "TJP1", "CDH1", "CLDN5", "CLDN11"))
# ggsave("umap_neuroepithelial_markers.png", height = 20, width = 30, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("NES", "VIM", "PAX6", "SOX2", "SOX10", "NOTCH1", "HES1", "OCLN", "TJP1", "CDH1", "CLDN5", "CLDN11"), 
        pt.size = 0)
# ggsave("violin_neuroepithelial_markers.png", height = 15, width = 35, dpi = 300)

# radial glia
FeaturePlot(seurat_dev, 
            features = c("VIM", "NES", "SOX2", "PAX6", "HES1", "HES5", "CDH2", "TNC", "GFAP", "FABP7", "SLC1A3", "SLC1A2"))
ggsave("umap_radialglia_markers.png", height = 20, width = 30, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("VIM", "NES", "SOX2", "PAX6", "HES1", "HES5", "CDH2", "TNC", "GFAP", "FABP7", "SLC1A3", "SLC1A2"), 
        pt.size = 0)
# ggsave("violin_radialglia_markers.png", height = 15, width = 35, dpi = 300)

# neural progenitor cell (NPC)
FeaturePlot(seurat_dev, 
            features = c("TBR1", "EOMES", "ASCL1", "DCX", "NEUROD1", "STMN1", "HES6", "DLL1", "DLL3"))
# ggsave("umap_npc_markers.png", height = 20, width = 25, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("TBR1", "EOMES", "ASCL1", "DCX", "NEUROD1", "STMN1", "HES6", "DLL1", "DLL3"), 
        pt.size = 0)
# ggsave("violin_npc_markers.png", height = 15, width = 35, dpi = 300)

# neuron
FeaturePlot(seurat_dev, 
            features = c("RBFOX3", "MAP2", "NEFL", "SYP", "GAP43", "TUBB3", "DLG4", "SYN1"))
ggsave("umap_neuron_markers.png", height = 20, width = 25, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("RBFOX3", "MAP2", "NEFL", "SYP", "GAP43", "TUBB3", "DLG4", "SYN1"), 
        pt.size = 0)
ggsave("violin_neuron_markers.png", height = 15, width = 30, dpi = 300)

# astrocyte
FeaturePlot(seurat_dev, 
            features = c("GFAP", "SLC1A2", "NEAT1", "GLUL", "S100B", "ALDH1L1", "AQP4", "GJA1"))
ggsave("umap_astrocyte_markers.png", height = 20, width = 25, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("GFAP", "SLC1A2", "NEAT1", "GLUL", "S100B", "ALDH1L1", "AQP4", "GJA1"), 
        pt.size = 0)
ggsave("violin_astrocyte_markers.png", height = 15, width = 30, dpi = 300)

# oligodendrocyte precusor cell (OPC)
FeaturePlot(seurat_dev, 
            features = c("PDGFRA", "CSPG4", "OLIG1", "OLIG2", "SOX10"))
ggsave("umap_OPC_markers.png", height = 20, width = 15, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("PDGFRA", "CSPG4", "OLIG1", "OLIG2", "SOX10"), 
        pt.size = 0)
ggsave("violin_OPC_markers.png", height = 10, width = 30, dpi = 300)

# oligodendrocyte
FeaturePlot(seurat_dev, 
            features = c("CLDN11", "MBP", "OLIG1", "OLIG2", "MOG", "CNP", "PLP1"))
ggsave("umap_oligo_markers.png", height = 20, width = 25, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("CLDN11", "MBP", "OLIG1", "OLIG2", "MOG", "CNP", "PLP1"), 
        pt.size = 0)
ggsave("violin_oligo_markers.png", height = 15, width = 30, dpi = 300)

# microglia
FeaturePlot(seurat_dev, 
            features = c("TMEM119", "ITGAM", "PTPRC", "AIF1", "CX3CR1", "P2RY12", "TREM2", "CSF1R"))
ggsave("umap_microglia_markers.png", height = 20, width = 25, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("TMEM119", "ITGAM", "PTPRC", "AIF1", "CX3CR1", "P2RY12", "TREM2", "CSF1R"), 
        pt.size = 0)
ggsave("violin_microglia_markers.png", height = 15, width = 30, dpi = 300)

# tanycyte
FeaturePlot(seurat_dev, 
            features = c("VIM", "NES", "SOX2", "RAX", "SLC1A3", "FGF10", "FGFRR3", "FABP7", "MAFB", "SCN7A", "DIO2", "FRZB"))
ggsave("umap_tanycyte_markers.png", height = 20, width = 30, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("VIM", "NES", "SOX2", "RAX", "SLC1A3", "FGF10", "FGFRR3", "FABP7", "MAFB", "SCN7A", "DIO2", "FRZB"), 
        pt.size = 0)
ggsave("violin_tanycyte_markers.png", height = 15, width = 35, dpi = 300)

# ependymal cell
FeaturePlot(seurat_dev, 
            features = c("FOXJ1", "TUBA1A", "DNAI1", "DNAH5", "PROM1"))
ggsave("umap_ependymal_markers.png", height = 20, width = 15, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("FOXJ1", "TUBA1A", "DNAI1", "DNAH5", "PROM1"), 
        pt.size = 0)
ggsave("violin_ependymal_markers.png", height = 10, width = 30, dpi = 300)

# endothelial
FeaturePlot(seurat_dev, 
            features = c("CDH5", "CLDN5", "PECAM1", "VWF", "CD34", "MCAM"))
ggsave("umap_endothelial_markers.png", height = 20, width = 15, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("CDH5", "CLDN5", "PECAM1", "VWF", "CD34", "MCAM"), 
        pt.size = 0)
ggsave("violin_endothelial_markers.png", height = 10, width = 30, dpi = 300)

# vascular leptomenegial cell (VLMC)
FeaturePlot(seurat_dev, 
            features = c("ACTA2", "MYH11", "TAGLN", "CNN1"))
ggsave("umap_vlmc_markers.png", height = 15, width = 15, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("PDGFRA", "COL1A1", "COL1A2", "DCN", "LUM"), 
        pt.size = 0)
ggsave("violin_vlmc_markers.png", height = 10, width = 30, dpi = 300)

# vascular smooth muscle cells (vsmcs) 
FeaturePlot(seurat_dev, 
            features = c("ACTA2", "MYH11", "TAGLN", "CNN1", "DES", "MCAM", "STMN"))
ggsave("umap_vsmc_markers.png", height = 20, width = 15, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("ACTA2", "MYH11", "TAGLN", "CNN1", "DES", "MCAM", "STMN"), 
        pt.size = 0)
ggsave("violin_vsmc_markers.png", height = 10, width = 30, dpi = 300)

# pericytes 
FeaturePlot(seurat_dev, 
            features = c("PDGFRB", "CSPG4", "RGS5", "ABCC9", "ANPEP"))
ggsave("umap_pericyte_markers.png", height = 20, width = 35, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("PDGFRB", "CSPG4", "RGS5", "ABCC9", "ANPEP"), 
        pt.size = 0)
ggsave("violin_pericyte_markers.png", height = 10, width = 30, dpi = 300)

# --------------------------------------------------------------------------------
# cell type identification using paper markers
# --------------------------------------------------------------------------------

# neuroepithelial
FeaturePlot(seurat_dev, 
            features = c("VIM", "HES1", "TJP1", "NES", "HMGA2", "ARHGAP28"))
ggsave("umap_neuroepithelial_paper.png", height = 20, width = 15, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("VIM", "HES1", "TJP1", "NES", "HMGA2", "ARHGAP28"), 
        pt.size = 0)
ggsave("violin_neuroepithelial_paper.png", height = 10, width = 30, dpi = 300)

# radial glia
FeaturePlot(seurat_dev, 
            features = c("VIM", "HMGA2", "COL2A1", "SOX2", "HES1", "FABP7", "SLC1A3", "TTYH1", "MKI67", "TOP2A", "CDK1"))
ggsave("umap_radialglia_paper.png", height = 20, width = 30, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("VIM", "HMGA2", "COL2A1", "SOX2", "HES1", "FABP7", "SLC1A3", "TTYH1", "MKI67", "TOP2A", "CDK1"), 
        pt.size = 0)
ggsave("violin_radialglia_paper.png", height = 15, width = 35, dpi = 300)

# neural progenitor cell (NPC)
FeaturePlot(seurat_dev, 
            features = c("ASCL1", "FABP7", "NES", "NBPF10", "NBPF15", "STMN2", "PLAGL1", "SFTA3", "EGFR", "OLIG1", "OLIG2",
                         "SYT1", "NEUROG2", "NEUROD1", "NHLH2"))
ggsave("umap_npc_paper.png", height = 30, width = 30, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("ASCL1", "FABP7", "NES", "NBPF10", "NBPF15", "STMN2", "PLAGL1", "SFTA3", "EGFR", "OLIG1", "OLIG2",
                     "SYT1", "NEUROG2", "NEUROD1", "NHLH2"), 
        pt.size = 0)
ggsave("violin_npc_paper.png", height = 25, width = 35, dpi = 300)

# neuron
FeaturePlot(seurat_dev, 
            features = c("STMN2", "DCX", "SYT1"))
ggsave("umap_neuron_paper.png", height = 20, width = 20, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("STMN2", "DCX", "SYT1"), 
        pt.size = 0)
ggsave("violin_neuron_paper.png", height = 5, width = 30, dpi = 300)

# astrocyte
FeaturePlot(seurat_dev, 
            features = c("FAM107A", "HOPX", "AQP4", "AGT", "GJA1", "GFAP", "ALDOC", "SLC1A2"))
ggsave("umap_astrocyte_paper.png", height = 20, width = 20, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("FAM107A", "HOPX", "AQP4", "AGT", "GJA1", "GFAP", "ALDOC", "SLC1A2"), 
        pt.size = 0)
ggsave("violin_astrocyte_paper.png", height = 15, width = 30, dpi = 300)

# oligodendrocyte precusor cell (OPC)
FeaturePlot(seurat_dev, 
            features = c("EGFR", "PDGFRA", "GPR75-ASB3", "APOD", "OLIG1", "OLIG2"))
ggsave("umap_OPC_paper.png", height = 25, width = 15, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("EGFR", "PDGFRA", "GPR75-ASB3", "APOD", "OLIG1", "OLIG2"), 
        pt.size = 0)
ggsave("violin_OPC_paper.png", height = 10, width = 30, dpi = 300)

# oligodendrocyte
FeaturePlot(seurat_dev, 
            features = c("BCAS1", "MBP", "MOG", "MAG", "PLP1"))
ggsave("umap_olig_paper.png", height = 20, width = 15, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("BCAS1", "MBP", "MOG", "MAG", "PLP1"), 
        pt.size = 0)
ggsave("violin_olig_paper.png", height = 10, width = 30, dpi = 300)

# microglia
FeaturePlot(seurat_dev, 
            features = c("AIF1", "CX3CR1", "C1QB"))
ggsave("umap_microglia_paper.png", height = 20, width = 20, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("AIF1", "CX3CR1", "C1QB"), 
        pt.size = 0)
ggsave("violin_microglia_paper.png", height = 5, width = 30, dpi = 300)

# tanycyte
FeaturePlot(seurat_dev, 
            features = c("CRYM", "RAX"))
ggsave("umap_tanycyte_paper.png", height = 10, width = 20, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("CRYM", "RAX"), 
        pt.size = 0)
ggsave("violin_tanycyte_paper.png", height = 5, width = 30, dpi = 300)

# ependymal cell
FeaturePlot(seurat_dev, 
            features = c("CCDC74A", "CCDC74B", "LRTOMT", "CCDC153", "FOXJ1"))
ggsave("umap_ependymal_paper.png", height = 20, width = 15, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("CCDC74A", "CCDC74B", "LRTOMT", "CCDC153", "FOXJ1"), 
        pt.size = 0)
ggsave("violin_ependymal_paper.png", height = 10, width = 30, dpi = 300)

# endothelial
FeaturePlot(seurat_dev, 
            features = c("CLDN5", "FLT1", "SLC38A5", "ITM2A"))
ggsave("umap_endothelial_paper.png", height = 20, width = 20, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("CLDN5", "FLT1", "SLC38A5", "ITM2A"), 
        pt.size = 0)
ggsave("violin_endothelial_paper.png", height = 10, width = 30, dpi = 300)

# vascular leptomenegial cell (VLMC)
FeaturePlot(seurat_dev, 
            features = c("PTGDS", "COL1A1", "LUM", "DCN"))
ggsave("umap_vlmc_paper.png", height = 20, width = 20, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("PTGDS", "COL1A1", "LUM", "DCN"), 
        pt.size = 0)
ggsave("violin_vlmc_paper.png", height = 10, width = 30, dpi = 300)

# vascular smooth muscle cells (vsmcs) 
FeaturePlot(seurat_dev, 
            features = c("ACTA2", "TAGLN"))
ggsave("umap_vsmc_paper.png", height = 20, width = 10, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("ACTA2", "TAGLN"), 
        pt.size = 0)
ggsave("violin_vsmc_paper.png", height = 5, width = 20, dpi = 300)

# pericytes 
FeaturePlot(seurat_dev, 
            features = c("KCNJ6", "ABCC9", "RGS5"))
ggsave("umap_pericyte_paper.png", height = 20, width = 20, dpi = 300)

VlnPlot(seurat_dev, 
        features = c("KCNJ6", "ABCC9", "RGS5"), 
        pt.size = 0)
ggsave("violin_pericyte_paper.png", height = 5, width = 30, dpi = 300)

# --------------------------------------------------------------------------------
# remaining clusters
# --------------------------------------------------------------------------------

# get list of cells in cluster 31 - ATP clusxter
remove <- WhichCells(seurat_dev, idents = c(20, 27, 31, 32, 39, 40, 41))

# remove cells 
seurat_dev <- subset(seurat_dev, cells = setdiff(Cells(seurat_dev), remove_cells))

# renaming cluster numbers to cell types
new_cluster_ids <- c("neuron", "astrocyte", "OPC", "astrocyte", "neuron", "microglia", "neuron", "NPC", "VLMC", "NE", "neuron", 
                     "tanycyte", "ependymal", "neuron", "olig", "neuron", "neuron", "neuron", "neuron", "mural", "astrocyte",
                     "endothelial", "neuron", "neuron", "RG (quiescent)", "neuron", "microglia", "microglia", "neuron", "RG (cycling)", "neuron",
                     "neuron", "RG (transitioning)", "microglia", "neuron", "mural", "neuron", "neuron", "neuron")

names(new_cluster_ids) <- levels(seurat_dev)
seurat_dev <- RenameIdents(seurat_dev, new_cluster_ids)

# visualising clusters after cell-type identification
seurat_dev@meta.data$"cell_type" <- as.factor(seurat_dev@active.ident)
DimPlot(seurat_dev, reduction = "umap", group.by="cell_type", label = TRUE) 

ggsave("umap_major_annotated.png", height = 7, width = 8)