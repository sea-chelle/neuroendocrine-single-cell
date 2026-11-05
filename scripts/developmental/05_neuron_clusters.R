# loading required packages
library(here)
library(Seurat)
library(harmony)
library(tidyverse)
library(scales)

# set working directory to figures
setwd(here("plots"))

# subsetting neuronal clusters for subclustering
neuron_dev <- subset(seurat_dev, idents = "Neuron")

# normalise and final variable features
neuron_dev <- NormalizeData(neuron_dev)
neuron_dev <- FindVariableFeatures(neuron_dev,selection.method = "vst", nfeatures = 2000)

LabelPoints(plot = VariableFeaturePlot(neuron_dev),
            points = head(VariableFeatures(neuron_dev), 10),
            repel = TRUE, xnudge = 0, ynudge = 0)

# continue with standard seurat workflow
neuron_dev <- ScaleData(neuron_dev, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
neuron_dev <- RunPCA(neuron_dev, features = VariableFeatures(object = neuron_dev))
ElbowPlot(neuron_dev, reduction = "pca", ndims = 50)

# run harmony for batch correction
neuron_dev <- RunHarmony(neuron_dev, group.by.vars = "dataset")

# UMAP dimensionality reduction clusters
neuron_dev <- FindNeighbors(neuron_dev, reduction = "harmony", dims = 1:40)
neuron_dev <- FindClusters(neuron_dev, resolution = 0.5)
neuron_dev <- RunUMAP(neuron_dev, reduction = "harmony", dims = 1:40)

# visualising clusters
DimPlot(neuron_dev, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
ggsave("UMAP_neuron_unannotated.png", height = 7, width = 8)

# looking at cluster sizes
table(neuron_dev@meta.data$seurat_clusters) 

# save unannotated seurat object
saveRDS(neuron_dev, here("r_obj", "neuron_dev_unannotated.rds"))

# glutamageric markers
FeaturePlot(neuron_dev, 
            features = c("SLC17A6", "SLC17A8", "GLS"))
ggsave("umap_glutamatergic_marker.png", height = 10, width = 10, dpi = 300)

VlnPlot(neuron_dev, 
        features = c("SLC17A6", "SLC17A8", "GLS"),
        pt.size = 0)
ggsave("violin_glutamatergic_marker.png", height = 5, width = 25, dpi = 300)

# GABAergic markers
FeaturePlot(neuron_dev, 
            features = c("SLC32A1", "SLC6A1", "GAD1", "GAD2", "DLX1", "DLX2", "DLX6", "ARX"))
ggsave("umap_gabaergic_marker.png", height = 10, width = 12, dpi = 300)

VlnPlot(neuron_dev, 
        features = c("SLC32A1", "SLC6A1", "GAD1", "GAD2", "DLX1", "DLX2", "DLX6", "ARX"),
        pt.size = 0)
ggsave("violin_gabaergic_marker.png", height = 10, width = 30, dpi = 300)

# histamergic 
FeaturePlot(neuron_dev, 
            features = "HDC")
ggsave("umap_histamatergic_marker.png", height = 5, width = 5)

VlnPlot(neuron_dev, 
        features = "HDC",
        pt.size = 0)
ggsave("violin_histamatergic_marker.png", height = 5, width = 10, dpi = 300)

# immature neurons
FeaturePlot(neuron_dev, 
            features = c("STMN2", "TUBB3", "DCX", "GAP43"))
ggsave("umap_immatureneuron_marker.png", height = 10, width = 10)

VlnPlot(neuron_dev, 
        features = c("STMN2", "TUBB3", "DCX", "GAP43"),
        pt.size = 0)
ggsave("violin_immatureneuron_marker.png", height = 8, width = 25, dpi = 300)

# mature neurons
FeaturePlot(neuron_dev, 
            features = c("MAP2", "SYT1", "SNAP25", "SYN1", "RBFOX3", "NEFL", "NEFM"))
ggsave("umap_matureneuron_marker.png", height = 12, width = 15)

VlnPlot(neuron_dev, 
        features = c("MAP2", "SYT1", "SNAP25", "SYN1", "RBFOX3", "NEFL", "NEFM"),
        pt.size = 0)
ggsave("violin_matureneuron_marker.png", height = 10, width = 25, dpi = 300)

# find differentially expressed genes from each cell type
neuron_dev_markers <- FindAllMarkers(neuron_dev, 
                                     only.pos = TRUE, 
                                     min.pct = 0.25,
                                     min.diff.pct = 0.25,
                                     logfc.threshold = 0.25)

# get top 10 markers for each cluster
top10_markers <- function(x){
  neuron_dev_markers[neuron_dev_markers$cluster == x, ] %>% head(n = 10)
  }

# create df of top 10 markers for each cluster
top10_markers_neuron_df <- map_dfr(0:32, top10_markers)

# save as .csv file
write.csv(top10_markers_neuron_df, here("csv", "top10_markers_neuron_dev.csv"))

# renaming cluster numbers to cell types
new_cluster_ids_neuron <- c("Immature", "Immature", "Glutamatergic", "Immature", "GABAergic", "GABAergic", "Glutamatergic", "GABAergic", "GABAergic", "GABAergic",
                            "Histaminergic", "GABAergic", "Glutamatergic", "GABAergic", "GABAergic", "Other", "GABAergic", "Glutamatergic", "Glutamatergic", "Glutamatergic",
                            "Other", "Other", "Histaminergic", "GABAergic", "Glutamatergic", "Other", "GABAergic", "Glutamatergic", "Glutamatergic", "GABAergic",
                            "Glutamatergic", "Other", "GABAergic")

names(new_cluster_ids_neuron) <- levels(neuron_dev)
neuron_dev <- RenameIdents(neuron_dev, new_cluster_ids_neuron)

# visualising clusters after cell-type identification
neuron_dev@meta.data$"cell_type" <- as.factor(neuron_dev@active.ident)
table(neuron_dev@meta.data$cell_type) # looking at cluster sizes

DimPlot(neuron_dev, reduction = "umap", group.by = "cell_type", label = FALSE, pt.size = 1.5) +
  scale_colour_manual(values = c("Glutamatergic" = "#78A8D1", "GABAergic" = "#A995C9", "Other" = "#D78FA9", 
                                 "Histaminergic" = "#5FA9A2", "Immature" = "grey80")) +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  theme_void() +
  theme(plot.title = element_blank(), 
        legend.text = element_text(size = 7),
        legend.key.height = unit(4, "mm"),
        legend.box.spacing = unit(0, "pt"))

ggsave("umap_neuron_annotated.tiff", height = 2.5, width = 3.3, dpi = 600, compression = "lzw")

# set order
neuron_dev$cell_type <- factor(neuron_dev$cell_type,
                               levels = c("Immature", "Glutamatergic", "GABAergic", "Histaminergic", "Other"))

# plot dotplot of markers
DotPlot(neuron_dev, assay = "RNA", group.by = "cell_type",
        features = c("HDC", # histaminergic
          "ARX", "DLX2","GAD2", "SLC32A1", # GABAergic
          "SLC17A6", # glutamatergic
          "SNAP25", "MAP2", "SYP", "STMN2", "TUBB3", # neuronal
          "VIM", "HES1" # immature neuronal
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
        panel.heights = unit(13 * 0.3, "cm"),
        panel.widths = unit(length(unique(neuron_dev$cell_type)) * 0.5, "cm"))

ggsave("major_marker_dotplot.tiff", height = 2.6, width = 2.1, dpi = 600, compression = "lzw")

# save annotated Seruat object
saveRDS(neuron_dev, here("r_obj", "neuron_dev.rds"))