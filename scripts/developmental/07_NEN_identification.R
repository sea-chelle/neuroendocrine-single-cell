# load required packaged
library(here)
library(tidyverse)
library(Seurat)

# set working directory
setwd(here("plots"))

# subset PVN clusters for subclustering
PVN_dev <- subset(nuclei_dev, idents = c("PVN"))

# normalise and scale data as usual
PVN_dev <- NormalizeData(PVN_dev)
PVN_dev <- FindVariableFeatures(PVN_dev, selection.method = "vst",  nfeatures = 2000)

# remove neuropeptides from highly variable genes
peptides <- c("CRH", "OXT", "AVP", "TRH", "SST", "GHRH", "POMC", "NPY", "AGRP", "KISS1", "TAC3", "GNRH1")
hvgs_pvn <- setdiff(VariableFeatures(PVN_dev), peptides)

# continue with standard seurat workflow
PVN_dev <- ScaleData(PVN_dev, features = hvgs_pvn, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
PVN_dev <- RunPCA(PVN_dev, features = hvgs_pvn)
ElbowPlot(PVN_dev)

# UMAP dimensionality reduction clusters
PVN_dev <- FindNeighbors(PVN_dev, dims = 1:20)
PVN_dev <- FindClusters(PVN_dev, resolution = 0.8)
PVN_dev <- RunUMAP(PVN_dev, dims = 1:20)

# visualising clusters
DimPlot(PVN_dev, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
ggsave("UMAP_PVN_unannotated.png", height = 5, width = 6)

# cluster sizes
table(PVN_dev@meta.data$seurat_clusters) 

# neuropeptidergic populations in the PVN
FeaturePlot(PVN_dev, features = c("CRH", "AVP", "OXT", "TRH"))
ggsave("UMAP_PVN_neuropeptides.png", height = 15, width = 15, dpi = 300)

VlnPlot(PVN_dev, features = c("CRH", "AVP", "OXT", "TRH"))
ggsave("vln_PVN_neuropeptides.png", height = 10, width = 30)

# neuroendocrine co-markers in the PVN
FeaturePlot(PVN_dev, features = c("NR3C1", "SCGN", "AGTR1", "TAC1", "TH", "ADARB2", "CALB1", "KCNMB4", "FOXP1",
                                  "THRB", "NFIB", "CALB2", "SFRP2", "OMP", "UCN3", "GHRH"))
ggsave("UMAP_PVN_neuroendocrine.png", height = 20, width = 25, dpi = 300)

VlnPlot(PVN_dev, features = c("NR3C1", "SCGN", "PCSK1", "PCSK2", "AGTR1", "TAC1", "TH", "ADARB2", "CALB1", "KCNMB4", "FOXP1",
                                   "THRB", "NFIB", "CALB2", "SFRP2","OMP", "UCN3", "GHRH"))
ggsave("vln_PVN_neuroendocrine.png", height = 20, width = 25)

# subset POA clusters for subclustering
POA_dev <- subset(nuclei_dev, idents = c("POA"))

# normalise and scale data as usual
POA_dev <- NormalizeData(POA_dev)
POA_dev <- FindVariableFeatures(POA_dev, selection.method = "vst",  nfeatures = 2000)

# remove neuropeptides from highly variable genes
hvgs_PO <- setdiff(VariableFeatures(POA_dev), peptides)

# continue with standard seurat workflow
POA_dev <- ScaleData(POA_dev, features = hvgs_PO, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
POA_dev <- RunPCA(POA_dev, features = hvgs_PO)
ElbowPlot(POA_dev)

# UMAP dimensionality reduction clusters
POA_dev <- FindNeighbors(POA_dev, dims = 1:20)
POA_dev <- FindClusters(POA_dev, resolution = 0.8)
POA_dev <- RunUMAP(POA_dev, dims = 1:20)

# visualise clusters
DimPlot(POA_dev, reduction = "umap", group.by = "seurat_clusters", label = T)
ggsave("UMAP_POA_unannotated.png", height = 5, width = 6)

# cluster sizes
table(POA_dev@meta.data$seurat_clusters) 

# identify neuropeptidergic populations in the POA
FeaturePlot(POA_dev, features = c("GNRH1"))
ggsave("UMAP_POA_neuropeptides.png", height = 10, width = 10, dpi = 300)

VlnPlot(POA_dev, features = c("GNRH1"))
ggsave("vln_POA_neuropeptides.png", height = 5, width = 10)

# neuroendocrine co-markers in the POA
FeaturePlot(POA_dev, features = c("KISS1R", "COCH", "CXCL14", "EBF3", "SEMA3C", "SEMA4D", "SOX4", 
                                  "PDE11A", "FREM1", "SPAG16", "MEIS1", "CFAP206", "KCNMB2", "EBF3", "MCTP2", "SIX6",
                                  "ISL1"))
ggsave("UMAP_POA_neuroendocrine.png", height = 25, width = 30, dpi = 300)

VlnPlot(POA_dev, features = c("KISS1R", "COCH", "CXCL14", "EBF3", "SEMA3C", "SEMA4D", "SOX4", 
                              "PDE11A", "FREM1", "SPAG16", "MEIS1", "CFAP206", "KCNMB2", "EBF3", "MCTP2", "SIX6",
                              "ISL1"))
ggsave("vln_POA_neuroendocrine.png", height = 15, width = 35)

# subset anterior hypothalamus clusters for subclustering
anterior_dev <- subset(nuclei_dev, idents = c("AH"))

# normalise and scale data as usual
anterior_dev <- NormalizeData(anterior_dev)
anterior_dev <- FindVariableFeatures(anterior_dev, selection.method = "vst",  nfeatures = 2000)

# remove neuropeptides from highly variable genes
hvgs_AH <- setdiff(VariableFeatures(anterior_dev), peptides)

# continue with standard seurat workflow
anterior_dev <- ScaleData(anterior_dev, features = hvgs_AH, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
anterior_dev <- RunPCA(anterior_dev, features = hvgs_AH)
ElbowPlot(anterior_dev)

# UMAP dimensionality reduction clusters
anterior_dev <- FindNeighbors(anterior_dev, dims = 1:20)
anterior_dev <- FindClusters(anterior_dev, resolution = 0.8)
anterior_dev <- RunUMAP(anterior_dev, dims = 1:20)

# visualise clusters
DimPlot(anterior_dev, reduction = "umap", group.by = "seurat_clusters", label = T)
ggsave("UMAP_AH_unannotated.png", height = 5, width = 6)

# cluster sizes
table(anterior_dev@meta.data$seurat_clusters) 

# neuropeptidergic populations in the AH
FeaturePlot(anterior_dev, features = c("CRH", "AVP", "OXT", "TRH", "GNRH1"))
ggsave("UMAP_AH_neuropeptides.png", height = 25, width = 15, dpi = 300)

VlnPlot(anterior_dev, features = c("CRH", "AVP", "OXT", "TRH", "GNRH1"))
ggsave("vln_AH_neuropeptides.png", height = 15, width = 30)

# neuroendocrine co-markers in the AH
FeaturePlot(anterior_dev, features = c("NR3C1", "SCGN", "AGTR1", "TAC1", "TH", "ADARB2", "CALB1", "KCNMB4", "FOXP1",
                                       "THRB", "NFIB", "CALB2", "SFRP2", "OMP", "UCN3", "GHRH"))
ggsave("UMAP_AH_PVN_neuroendocrine.png", height = 20, width = 25, dpi = 300)

VlnPlot(anterior_dev, features = c("NR3C1", "SCGN", "PCSK1", "PCSK2", "AGTR1", "TAC1", "TH", "ADARB2", "CALB1", "KCNMB4", "FOXP1",
                                   "THRB", "NFIB", "CALB2", "SFRP2","OMP", "UCN3", "GHRH"))
ggsave("vln_AH_PVN_neuroendocrine.png", height = 20, width = 25)

FeaturePlot(anterior_dev, features = c("KISS1R", "COCH", "CXCL14", "EBF3", "SEMA3C", "SEMA4D", "SOX4", 
                                       "PDE11A", "FREM1", "SPAG16", "MEIS1", "CFAP206", "KCNMB2", "EBF3", "MCTP2", "SIX6",
                                       "ISL1"))
ggsave("UMAP_AH_PO_neuroendocrine.png", height = 25, width = 30, dpi = 300)

VlnPlot(anterior_dev, features = c("KISS1R", "COCH", "CXCL14", "EBF3", "SEMA3C", "SEMA4D", "SOX4", 
                                   "PDE11A", "FREM1", "SPAG16", "MEIS1", "CFAP206", "KCNMB2", "EBF3", "MCTP2", "SIX6",
                                   "ISL1"))
ggsave("vln_AH_PO_neuroendocrine.png", height = 15, width = 35)

# subset ARC clusters for subclustering
ARC_dev <- subset(nuclei_dev, idents = c("ARC"))

# normalise and scale data as usual
ARC_dev <- NormalizeData(ARC_dev)
ARC_dev <- FindVariableFeatures(ARC_dev,  selection.method = "vst",  nfeatures = 2000)

# remove neuropeptides from highly variable genes
hvgs_ARC <- setdiff(VariableFeatures(ARC_dev), peptides)

# continue with standard seurat workflow
ARC_dev <- ScaleData(ARC_dev, features = hvgs_ARC, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
ARC_dev <- RunPCA(ARC_dev, features = hvgs_ARC)
ElbowPlot(ARC_dev)

# UMAP dimensionality reduction clusters
ARC_dev <- FindNeighbors(ARC_dev, dims = 1:20)
ARC_dev <- FindClusters(ARC_dev, resolution = 0.8)
ARC_dev <- RunUMAP(ARC_dev, dims = 1:20)

# visualise clusters
DimPlot(ARC_dev, reduction = "umap", group.by = "seurat_clusters", label = T)

# cluster sizes
table(ARC_dev@meta.data$seurat_clusters) 
ggsave("UMAP_ARC_unannotated.png", height = 5, width = 6)

# neuropeptidergic populations in the ARC
FeaturePlot(ARC_dev, features = c("GHRH", "TH"))
ggsave("UMAP_ARC_neuropeptides.png", height = 8, width = 15, dpi = 300)

VlnPlot(ARC_dev, features = c("GHRH", "TH"))
ggsave("vln_ARC_neuropeptides.png", height = 5, width = 15)

# neuroendocrine co-markers in the ARC
FeaturePlot(ARC_dev, features = c("TRH", "GAL",
                                  "SLC6A3", "DDC", "SLC18A2", "GSX1", "FOXP2", "DLX1", "DLX2",
                                  "PRLR", "ARX", "LHX1", "SST"))
ggsave("UMAP_ARC_neuroendocrine.png", height = 25, width = 30, dpi = 300)

VlnPlot(ARC_dev, features = c("TRH", "GAL",
                              "SLC6A3", "DDC", "SLC18A2", "GSX1", "FOXP2", "DLX1", "DLX2",
                              "PRLR", "ARX", "LHX1", "SST"))
ggsave("vln_ARC_neuroendocrine.png", height = 15, width = 35)

# subsetting tuberal hypothalamus clusters for subclustering
tuberal_dev <- subset(nuclei_dev, idents = c("Tuberal"))

# normalise and scale data as usual
tuberal_dev <- NormalizeData(tuberal_dev)
tuberal_dev <- FindVariableFeatures(tuberal_dev, selection.method = "vst",  nfeatures = 2000)
hvgs_unknown <- setdiff(VariableFeatures(tuberal_dev), peptides)

# continue with standard seurat workflow
tuberal_dev <- ScaleData(tuberal_dev, features = hvgs_unknown, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
tuberal_dev <- RunPCA(tuberal_dev, features = hvgs_unknown)
ElbowPlot(tuberal_dev)

# UMAP dimensionality reduction clusters
tuberal_dev <- FindNeighbors(tuberal_dev, dims = 1:20)
tuberal_dev <- FindClusters(tuberal_dev, resolution = 0.8)
tuberal_dev <- RunUMAP(tuberal_dev, dims = 1:20)

# visualise clusters
DimPlot(tuberal_dev, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
ggsave("UMAP_tuberal_unannotated.png", height = 5, width = 6)

# cluster sizes
table(tuberal_dev@meta.data$seurat_clusters) 

# neuropeptidergic populations in the tuberal hypothalamus
FeaturePlot(tuberal_dev, features = c("GHRH", "TH"))
ggsave("UMAP_tuberal_neuropeptides.png", height = 8, width = 15, dpi = 300)

VlnPlot(tuberal_dev, features = c("GHRH", "TH"))
ggsave("Vln_tuberal_neuropeptides.png", height = 5, width = 15)

# neuroendocrine co-markers in the tuberal hypothalamus
FeaturePlot(tuberal_dev, features = c("TRH", "GAL", "SLC6A3", "DDC", "SLC18A2", "GSX1", "FOXP2", "DLX1", "DLX2", "PRLR", "ARX", "LHX1", "SST"))
ggsave("UMAP_tuberal_neuroendocrine.png", height = 25, width = 30, dpi = 300)

VlnPlot(tuberal_dev, features = c("TRH", "GAL", "SLC6A3", "DDC", "SLC18A2", "GSX1", "FOXP2", "DLX1", "DLX2", "PRLR", "ARX", "LHX1", "SST"))
ggsave("vln_tuberal_neuroendocrine.png", height = 15, width = 35)

# check number and proportion of rare populations in each cluster
table(PVN_dev$seurat_clusters[WhichCells(PVN_dev, expression = CRH > 0, slot = "counts")])
100 * table(PVN_dev$seurat_clusters[WhichCells(PVN_dev, expression = CRH > 0, slot = "counts")]) /
  table(PVN_dev$seurat_clusters) # CRH in the PVN

table(POA_dev$seurat_clusters[WhichCells(POA_dev, expression = GNRH1 > 0, slot = "counts")])
100 * table(POA_dev$seurat_clusters[WhichCells(POA_dev, expression = GNRH1 > 0, slot = "counts")]) /
  table(POA_dev$seurat_clusters) # GNRH1 in the POA

table(ARC_dev$seurat_clusters[WhichCells(ARC_dev, expression = DDC > 0, slot = "counts")])
100 * table(ARC_dev$seurat_clusters[WhichCells(ARC_dev, expression = DDC > 0, slot = "counts")]) /
  table(ARC_dev$seurat_clusters) # TIDA neurons in the ARC

# create NEN vs non_NEN column in metadata
nuclei_dev$NEN_type <- "non_NEN" # label every cell as non_NEN before rewriting NEN labels

# add labels for NENs
nuclei_dev$NEN_type[colnames(nuclei_dev) %in% 
                      colnames(PVN_dev)[PVN_dev$seurat_clusters %in% c(5, 8, 10, 13)]] <- "NEN"

nuclei_dev$NEN_type[colnames(nuclei_dev) %in% 
                      colnames(ARC_dev)[ARC_dev$seurat_clusters %in% c(6, 12)]] <- "NEN"

nuclei_dev$NEN_type[colnames(nuclei_dev) %in% 
                      colnames(tuberal_dev)[tuberal_dev$seurat_clusters %in% c(8)]] <- "NEN"


# number of NENs there are
table(nuclei_dev$NEN_type)

# add labels for NEN subtypes
nuclei_dev$NEN_subtype <- "non_NEN" 

nuclei_dev$NEN_subtype[colnames(nuclei_dev) %in% 
                         colnames(PVN_dev)[PVN_dev$seurat_clusters %in% 5]] <- "AVP-CRH"

nuclei_dev$NEN_subtype[colnames(nuclei_dev) %in% 
                         colnames(PVN_dev)[PVN_dev$seurat_clusters %in% 13]] <- "AVP-OXT"

nuclei_dev$NEN_subtype[colnames(nuclei_dev) %in% 
                         colnames(PVN_dev)[PVN_dev$seurat_clusters %in% c(8, 10)]] <- "TRH"

nuclei_dev$NEN_subtype[colnames(nuclei_dev) %in% 
                         colnames(ARC_dev)[ARC_dev$seurat_clusters %in% c(6, 12)]] <- "GHRH"

nuclei_dev$NEN_subtype[colnames(nuclei_dev) %in% 
                         colnames(tuberal_dev)[tuberal_dev$seurat_clusters %in% 8]] <- "GHRH"

# number of each NEN subtype
table(nuclei_dev$NEN_subtype)

# pie chart for global proportion of NENs and non-NENs within hypothalamic nuclei
NEN_summary <- nuclei_dev@meta.data %>%
  count(NEN_type, name = "cell_count")

NEN_summary_pie <- NEN_summary %>%
  mutate(percentage = cell_count / sum(cell_count) * 100)

ggplot(NEN_summary_pie, aes(x = "", y = percentage, fill = NEN_type)) +
  geom_col(width = 1) +
  geom_text(aes(x = 1.6, label = paste0(round(percentage, 1), "%")),
            position = position_stack(vjust = 0.5),
            size = 2.5
  ) +
  coord_polar(theta = "y", clip = "off") +
  scale_fill_manual(values = c("NEN" = "#A995C9", "non_NEN" = "grey85")) +
  theme_void() +
  theme(legend.position = "none")

ggsave("NEN_pie_proportions.tiff", height = 1.5, width = 1.5, dpi = 600, compression = "lzw")

# stacked bar plot per region of NENs vs non-NENs
NEN_region_summary <- nuclei_dev@meta.data %>%
  filter(cell_type %in% c("PVN", "ARC", "Tuberal")) %>%
  count(cell_type, NEN_type, name = "cell_count") %>%
  group_by(cell_type) %>%
  mutate(percentage = 100 * cell_count / sum(cell_count)) %>%
  ungroup()

NEN_region_summary$cell_type <- factor(NEN_region_summary$cell_type, levels = c("PVN", "ARC", "Tuberal"))

ggplot(NEN_region_summary, aes(x = cell_type, y = percentage, fill = NEN_type)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_fill_manual( values = c("NEN" = "#A995C9", "non_NEN" = "grey85")) +
  labs(x = NULL,y = "Percentage of cells (%)") +
  theme_classic() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.line = element_line(linewidth = 0.3),
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        legend.box.spacing = unit(0, "pt"),
        legend.key.size = unit(0.3, "cm"))

ggsave("NEN_per_region.tiff", height = 1.5, width = 2.6, dpi = 600, compression = "lzw")

# set order
nuclei_dev$NEN_subtype <- factor(nuclei_dev$NEN_subtype,
                               levels = c("AVP-CRH", "AVP-OXT", "TRH", "GHRH", "non_NEN"))

# violin plot of neuropeptides and markers
NEN_violin <- VlnPlot(subset(nuclei_dev, subset = NEN_type == "NEN"),
  features = c("AVP", "CRH", "OXT", "TRH", "GHRH",
               "FOXP1", "NR3C1", "ADARB2", "KCNMB4", "NFIB", "GAL", "FOXP2", "GSX1"),
  group.by = "NEN_subtype"
  cols = c(rep("#A995C9", 5), rep("#5FA9A2", 8)),
  pt.size = 0, stack = TRUE, flip = TRUE) +
  NoLegend() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.title.x = element_blank(),
        panel.border = element_blank(),
        strip.text.y.right = element_text(size = 8, face = "plain"))

NEN_violin$layers[[1]]$aes_params$linewidth <- 0.3
NEN_violin

ggsave("NEN_violin.tiff", height = 4, width = 3, dpi = 600, compression = "lzw")
