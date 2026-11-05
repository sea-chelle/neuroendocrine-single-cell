# load required packaged
library(here)
library(tidyverse)
library(Seurat)
library(patchwork)
library(ggrepel)

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

# exlcuding preautonomic ADARB2+ neuropeptidergic clusters
VlnPlot(PVN_dev, features = c("ADARB2", "RELN", "SNCA", "NTNG1", "CACNA1G"),
  group.by = "seurat_clusters", pt.size = 0, stack = TRUE, flip = TRUE)

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

# subsetting unresolved hypothalamus clusters for subclustering
unresolved_dev <- subset(nuclei_dev, idents = c("Unresolved"))

# normalise and scale data as usual
unresolved_dev <- NormalizeData(unresolved_dev)
unresolved_dev <- FindVariableFeatures(unresolved_dev, selection.method = "vst",  nfeatures = 2000)
hvgs_unknown <- setdiff(VariableFeatures(unresolved_dev), peptides)

# continue with standard seurat workflow
unresolved_dev <- ScaleData(unresolved_dev, features = hvgs_unknown, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
unresolved_dev <- RunPCA(unresolved_dev, features = hvgs_unknown)
ElbowPlot(unresolved_dev)

# UMAP dimensionality reduction clusters
unresolved_dev <- FindNeighbors(unresolved_dev, dims = 1:20)
unresolved_dev <- FindClusters(unresolved_dev, resolution = 0.8)
unresolved_dev <- RunUMAP(unresolved_dev, dims = 1:20)

# visualise clusters
DimPlot(unresolved_dev, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
ggsave("UMAP_unresolved_unannotated.png", height = 5, width = 6)

# neuropeptidergic populations in the tuberal hypothalamus
FeaturePlot(unresolved_dev, features = c("AVP", "CRH", "OXT", "TRH", "GNRH1", "GHRH"))
ggsave("UMAP_unresolved_neuropeptides.png", height = 15, width = 15, dpi = 300)

VlnPlot(unresolved_dev, features = c("AVP", "CRH", "OXT", "TRH", "GNRH1", "GHRH"))
ggsave("Vln_unresolved_neuropeptides.png", height = 8, width = 15) # enrichment for PVN NENs - check lineage and co-markers

# check expression of PVN lineage markers
FeaturePlot(unresolved_dev, features = c("SIM1", "OTP", "POU3F2", "ARNT2"))
ggsave("UMAP_unresolved_neuropeptides.png", height = 15, width = 15, dpi = 300)

VlnPlot(unresolved_dev, features = c("SIM1", "OTP", "POU3F2", "ARNT2"))
ggsave("Vln_unresolved_neuropeptides.png", height = 8, width = 15) # enrichment in Clusters 2, 6, 13 and 15

# plot SIM1 and OTP on original UMAP to show enrichment in unknown clusters
FeaturePlot(nuclei_dev, features = c("SIM1", "OTP"), cols = c("grey80", "#7656A8"), keep.scale = "all", order = TRUE,  pt.size = 0.01) +
  plot_layout(guides = "collect") &
  guides(colour = guide_colourbar(title = "Expression")) &
  theme_void() &
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.height = unit(4, "mm"),
        legend.key.width = unit(4, "mm"), 
        legend.box.spacing = unit(0, "pt"))

ggsave("umap_SIM_OTP.tiff", height = 2.8, width = 4.5, dpi = 600, compression = "lzw")

# set idents as seurat clusters
unresolved_dev <- SetIdent(unresolved_dev, value = "seurat_clusters")

# calculate the median UMAP position of each cluster
cluster_centres <- as.data.frame(Embeddings(unresolved_dev, reduction = "umap"))

colnames(cluster_centres)[1:2] <- c("UMAP_1", "UMAP_2")

cluster_centres$cluster <- as.character(unresolved_dev$seurat_clusters)

cluster_centres <- cluster_centres %>%
  group_by(cluster) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

# create base FeaturePlots
plots <- FeaturePlot(unresolved_dev, features = c("SIM1", "OTP"), order = TRUE, pt.size = 0.01, combine = FALSE)

# add labels, identical colour scales and formatting 
plots <- lapply(plots, function(p) {
  p + geom_text_repel(data = cluster_centres,
      aes(x = UMAP_1, y = UMAP_2, label = cluster), inherit.aes = FALSE, size = 2.5, colour = "black", box.padding = 0.2, min.segment.length = 0) +
    scale_colour_gradient( low = "grey90", high = "#7656A8", limits = c(0, 4), oob = scales::squish, name = "Expression\nlevel") +
    theme_void() +
    theme(plot.title = element_text(size = 8, face = "bold",hjust = 0.5),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.key.height = unit(4, "mm"),
          legend.key.width = unit(4, "mm"),
          legend.box.spacing = unit(0, "pt"))
  })

# combine and collect the legends
unresolved_lineage_plot <- wrap_plots(plots, ncol = 2, guides = "collect") & 
  theme(legend.position = "right")

unresolved_lineage_plot

ggsave("unresolved_SIM1OTP_umap.tiff", unresolved_lineage_plot, height = 3.5, width = 6, dpi = 600, compression = "lzw")

# add PVN labels to these clusters in a new column in metadata
nuclei_dev$NEN_region <- as.character(nuclei_dev$cell_type)
PVN_cells <- Cells(unresolved_dev)[as.character(unresolved_dev$seurat_clusters) %in% c("2", "5", "6", "13", "15")]
nuclei_dev@meta.data[PVN_cells, "NEN_region"] <- "PVN"
nuclei_dev$NEN_region <- factor(nuclei_dev$NEN_region)

# neuroendocrine co-markers in the PVN
FeaturePlot(unresolved_dev, features = c("NR3C1", "SCGN", "AGTR1", "TAC1", "TH", "ADARB2", "CALB1", "KCNMB4", "FOXP1",
                                  "THRB", "NFIB", "CALB2", "SFRP2", "OMP", "UCN3", "GHRH"))
ggsave("UMAP_unresolved_neuroendocrine.png", height = 20, width = 25, dpi = 300)

VlnPlot(unresolved_dev, features = c("NR3C1", "SCGN", "PCSK1", "PCSK2", "AGTR1", "TAC1", "TH", "ADARB2", "CALB1", "KCNMB4", "FOXP1",
                              "THRB", "NFIB", "CALB2", "SFRP2","OMP", "UCN3", "GHRH"))
ggsave("vln_unresolved_neuroendocrine.png", height = 20, width = 25)

# exlcuding preautonomic ADARB2+ neuropeptidergic clusters
VlnPlot(unresolved_dev, features = c("ADARB2", "RELN", "SNCA", "NTNG1", "CACNA1G"),
        group.by = "seurat_clusters", pt.size = 0, stack = TRUE, flip = TRUE)

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

nuclei_dev$NEN_type[colnames(nuclei_dev) %in% 
                      colnames(unresolved_dev)[unresolved_dev$seurat_clusters %in% c(2, 6, 15)]] <- "NEN"

# number of NENs there are
table(nuclei_dev$NEN_type)

# add labels for NEN subtypes
nuclei_dev$NEN_subtype <- "non_NEN" 

nuclei_dev$NEN_subtype[colnames(nuclei_dev) %in% 
                         colnames(PVN_dev)[PVN_dev$seurat_clusters %in% 5]] <- "AVP-CRH"

nuclei_dev$NEN_subtype[colnames(nuclei_dev) %in% 
                         colnames(unresolved_dev)[unresolved_dev$seurat_clusters %in% 15]] <- "AVP-CRH"

nuclei_dev$NEN_subtype[colnames(nuclei_dev) %in% 
                         colnames(PVN_dev)[PVN_dev$seurat_clusters %in% 13]] <- "AVP-OXT"

nuclei_dev$NEN_subtype[colnames(nuclei_dev) %in% 
                         colnames(unresolved_dev)[unresolved_dev$seurat_clusters %in% 6]] <- "AVP-OXT"

nuclei_dev$NEN_subtype[colnames(nuclei_dev) %in% 
                         colnames(unresolved_dev)[unresolved_dev$seurat_clusters %in% 2]] <- "OXT"

nuclei_dev$NEN_subtype[colnames(nuclei_dev) %in% 
                         colnames(PVN_dev)[PVN_dev$seurat_clusters %in% c(8, 10)]] <- "TRH"

nuclei_dev$NEN_subtype[colnames(nuclei_dev) %in% 
                         colnames(ARC_dev)[ARC_dev$seurat_clusters %in% c(6, 12)]] <- "GHRH"

nuclei_dev$NEN_subtype[colnames(nuclei_dev) %in% 
                         colnames(tuberal_dev)[tuberal_dev$seurat_clusters %in% 8]] <- "GHRH"

# number of each NEN subtype
table(nuclei_dev$NEN_subtype)

# pie chart for global proportion of NENs and non-NENs within hypothalamic nuclei
NEN_subtype_summary <- nuclei_dev@meta.data %>%
  filter(NEN_type == "NEN", !is.na(NEN_subtype), NEN_subtype != "non_NEN") %>%
  count(NEN_subtype, name = "cell_count") %>%
  mutate(percentage = cell_count / sum(cell_count) * 100,
         NEN_subtype = factor(NEN_subtype, levels = c("AVP-CRH", "AVP-OXT", "OXT", "TRH", "GHRH"))) %>%
  arrange(NEN_subtype)

ggplot(NEN_subtype_summary, aes(x = "", y = percentage, fill = NEN_subtype)) +
  geom_col(width = 1) +
  geom_text(aes(x = 1.7, label = paste0(round(percentage, 1), "%")), position = position_stack(vjust = 0.5), size = 2.5) +
  coord_polar(theta = "y", clip = "off") +
  scale_fill_manual(values = c("AVP-CRH" = "#A995C9", "AVP-OXT" = "#78A8D1", "OXT" = "#D78FA9", "TRH" = "#C5A04A", "GHRH" = "#5FA9A2")) +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0, "pt"))

ggsave("NEN_pie_proportions.tiff", height = 2.5, width = 2.5, dpi = 600, compression = "lzw")

# set order
nuclei_dev$NEN_subtype <- factor(nuclei_dev$NEN_subtype,
                               levels = c("AVP-CRH", "AVP-OXT", "OXT", "TRH", "GHRH", "non_NEN"))

# violin plot of neuropeptides and markers
DotPlot(nuclei_dev, assay = "RNA", group.by = "NEN_subtype",
        features = c("GSX1", "FOXP2", "GAL", "NFIB", "KCNMB4", "CALB1", "TAC1", "FOXP1", "NR3C1", "GHRH", "TRH", "OXT", "CRH", "AVP"
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
        panel.heights = unit(10 * 0.5, "cm"),
        panel.widths = unit(length(unique(seurat_dev$cell_type)) * 0.3, "cm"))

ggsave("NEN_dotplot.tiff", height = 3.4, width = 2.9, dpi = 600, compression = "lzw")