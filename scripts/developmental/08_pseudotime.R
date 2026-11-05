# make cluster identity active for major cell types
seurat_dev <- SetIdent(seurat_dev, value = "seurat_clusters")

# progenitor clusters
progenitor_clusters <- c("6", "9", "10", "11", "13", "19", "21", "22", "25", "32", "39", "40")

# subset progenitor populations
progenitor_dev <- subset(seurat_dev, idents = progenitor_clusters)

# retain desired cluster order
progenitor_dev$seurat_clusters <- factor(progenitor_dev$seurat_clusters, levels = progenitor_clusters)

# specify marker sets for progenitors
trajectory_marker_sets <- list("Neuroepithelial/RG" = c("SOX2", "NES", "VIM", "HES1", "HES5", "FABP7", "SLC1A3", "TTYH1", "CDH2", "TJP1"),
                               "Cycling" = c("MKI67", "TOP2A", "PCNA", "CENPF", "UBE2C", "CKS2", "CKS1B"),
                               "Neurogenic" = c("ASCL1", "HES6", "DLL1", "NEUROG2", "NEUROD1", "NHLH1", "NHLH2", "DCX", "TUBB3"),
                               "Gliogenic/oligodendroglial" = c("SOX9", "NFIA", "NFIB", "EGFR", "OLIG1", "OLIG2", "PDGFRA", "CSPG4", "SOX10",
                                                                "S100B", "AQP4", "SPARCL1"))

# retain only genes present in the object
trajectory_marker_sets <- lapply(trajectory_marker_sets, function(x) intersect(x, rownames(progenitor_dev)))

# dot plot of progenitor markers within all progenitor clusters
DotPlot(progenitor_dev, features = trajectory_marker_sets, group.by = "seurat_clusters", assay = "RNA", scale = TRUE) +
  RotatedAxis() +
  labs(x = NULL, y = "Cluster", colour = "Average\nexpression", size = "Percent\nexpressed") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        strip.text = element_text(size = 7))

# include clusters 22 (neuroepithelial root), 25, 32, 39, 40 (cycling neural progenitors without gliogenic bias), 13 (IPCs without gliogenic bias)
# core neural progenitor clusters
progenitor_clusters <- c("22", "25", "32", "39", "40", "13")
progenitor_cells <- colnames(seurat_dev)[as.character(seurat_dev$seurat_clusters) %in% progenitor_clusters]

# set seed for reproducibility
set.seed(123)

# get metadata from immune cells
immature_meta <- immature_dev@meta.data %>%
  rownames_to_column("cell")

# sample up to 500 immature cells per sample × immature cluster (many more immature neurons than progenitors or mature neurons)
immature_sampled_cells <- immature_meta %>%
  group_by(dataset, sample, seurat_clusters) %>%
  group_modify(~ slice_sample(.x, n = min(500, nrow(.x)))) %>%
  ungroup() %>%
  pull(cell)

# see how many immature neurons were retained
length(immature_sampled_cells)

# see how many immature neurons were retained per dataset
table(immature_dev$dataset[colnames(immature_dev) %in% immature_sampled_cells])

# mature PVN neurons, including NEN and non-NEN PVN populations
pvn_cells <- colnames(nuclei_dev)[as.character(nuclei_dev$NEN_region) == "PVN"]

# combine unique cell barcodes
trajectory_cells <- unique(c(progenitor_cells, immature_sampled_cells, pvn_cells))

# ensure every barcode exists in the original object
trajectory_cells <- intersect(trajectory_cells, colnames(seurat_dev))

# subset seurat_dev for cells to include in PVN trajectory
PVN_trajectory_dev <- subset(seurat_dev, cells = trajectory_cells)

# add trajectory stage for each cluster
PVN_trajectory_dev$trajectory_stage <- NA_character_
PVN_trajectory_dev$trajectory_stage[as.character(PVN_trajectory_dev$seurat_clusters) == "22"] <- "Neuroepithelial/RG"
PVN_trajectory_dev$trajectory_stage[as.character(PVN_trajectory_dev$seurat_clusters) %in% c("25", "32", "39", "40")] <- "Cycling progenitor"
PVN_trajectory_dev$trajectory_stage[as.character(PVN_trajectory_dev$seurat_clusters) == "13"] <- "IPC"
PVN_trajectory_dev$trajectory_stage[colnames(PVN_trajectory_dev) %in% immature_sampled_cells] <- "Immature neuron"
PVN_trajectory_dev$trajectory_stage[colnames(PVN_trajectory_dev) %in% pvn_cells] <- "PVN neuron"

# set order
PVN_trajectory_dev$trajectory_stage <- factor(PVN_trajectory_dev$trajectory_stage,
                                              levels = c("Neuroepithelial/RG", "Cycling progenitor", "IPC", "Immature neuron", "PVN neuron"))

# see how many of each cell type there are
table(PVN_trajectory_dev$trajectory_stage, useNA = "ifany")

# see how many of each trajectory stage there are per GW
table(PVN_trajectory_dev$trajectory_stage, PVN_trajectory_dev$gestational.week)