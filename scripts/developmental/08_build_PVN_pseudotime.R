# load required packages
library(tidyverse)
library(here)
library(Seurat)
library(monocle3)
library(harmony)
library(SeuratWrappers)
library(patchwork)

set.seed(123)

# set Idents to seurat cluster number
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

# core neural progenitor clusters
progenitor_clusters <- c("22", "25", "32", "39", "40", "13")
progenitor_cells <- colnames(seurat_dev)[as.character(seurat_dev$seurat_clusters) %in% progenitor_clusters]

# sample up to 500 immature neurons per dataset, sample and neuronal cluster.
immature_dev <- subset(neuron_dev, subset = cell_type == "Immature")

immature_meta <- immature_dev@meta.data %>%
  rownames_to_column("cell") %>%
  arrange(dataset, sample, seurat_clusters, cell)

set.seed(123)
immature_sampled_cells <- immature_meta %>%
  group_by(dataset, sample, seurat_clusters) %>%
  group_modify(~ slice_sample(.x, n = min(500, nrow(.x)))) %>%
  ungroup() %>%
  pull(cell)

writeLines(sort(immature_sampled_cells),
           file.path(table_dir, "immature_sampled_cells.txt"))

# number of immature neurons retained
length(immature_sampled_cells)

# number of immature neurons retained per dataset
table(immature_dev$dataset[colnames(immature_dev) %in% immature_sampled_cells])

# retain all mature PVN cells from the reannotated nuclei object
pvn_cells <- colnames(nuclei_dev)[which(as.character(nuclei_dev$NEN_region) == "PVN")]

# subset for all cells to be included in trajectory
trajectory_cells <- unique(c(progenitor_cells, immature_sampled_cells, pvn_cells))
trajectory_cells <- intersect(trajectory_cells, colnames(seurat_dev))
PVN_trajectory_dev <- subset(seurat_dev, cells = trajectory_cells)

# add cell type category
PVN_trajectory_dev$trajectory_stage <- case_when(colnames(PVN_trajectory_dev) %in% pvn_cells ~ "PVN neuron",
                                                 colnames(PVN_trajectory_dev) %in% immature_sampled_cells ~ "Immature neuron",
                                                 as.character(PVN_trajectory_dev$seurat_clusters) == "22" ~ "Neuroepithelial/RG",
                                                 as.character(PVN_trajectory_dev$seurat_clusters) %in% c("25", "32", "39", "40") ~ "Cycling progenitor",
                                                 as.character(PVN_trajectory_dev$seurat_clusters) == "13" ~ "IPC", TRUE ~ NA_character_)

# set order
stage_order <- c("Neuroepithelial/RG", "Cycling progenitor", "IPC", "Immature neuron", "PVN neuron")
week_order <- c("GW6", "GW7", "GW8", "GW10", "GW12", "GW15", "GW16", "GW18", "GW19", "GW20", "GW22", "GW25")
PVN_trajectory_dev$trajectory_stage <- factor(PVN_trajectory_dev$trajectory_stage, levels = stage_order)
PVN_trajectory_dev$gestational.week <- factor(PVN_trajectory_dev$gestational.week, levels = week_order)

# transfer the nuclei annotations (progenitors and immature neurons remain NA)
matched_nuclei <- match(colnames(PVN_trajectory_dev), colnames(nuclei_dev))
PVN_trajectory_dev$NEN_region <- as.character(nuclei_dev$NEN_region[matched_nuclei])
PVN_trajectory_dev$NEN_type <- as.character(nuclei_dev$NEN_type[matched_nuclei])
PVN_trajectory_dev$NEN_subtype <- as.character(nuclei_dev$NEN_subtype[matched_nuclei])

PVN_trajectory_dev$trajectory_identity <- as.character(PVN_trajectory_dev$trajectory_stage)
mature_pvn <- PVN_trajectory_dev$trajectory_stage == "PVN neuron"
PVN_trajectory_dev$trajectory_identity[mature_pvn] <- paste0("PVN: ", PVN_trajectory_dev$NEN_subtype[mature_pvn])

# ensure DefaulAssay is RNA
DefaultAssay(PVN_trajectory_dev) <- "RNA"

# find variable features, rescale and rerun PCA
PVN_trajectory_dev <- FindVariableFeatures(PVN_trajectory_dev, selection.method = "vst", nfeatures = 2000)
PVN_trajectory_dev <- ScaleData(PVN_trajectory_dev, features = VariableFeatures(PVN_trajectory_dev), verbose = FALSE)
PVN_trajectory_dev <- RunPCA(PVN_trajectory_dev, features = VariableFeatures(PVN_trajectory_dev), npcs = 50, verbose = FALSE)

# rerun Harmony for batch correction
PVN_trajectory_dev <- RunHarmony(PVN_trajectory_dev, group.by.vars = "dataset")

# rerun UMAP for cells in trajectory
PVN_trajectory_dev <- RunUMAP(PVN_trajectory_dev, reduction = "harmony", dims = 1:20, seed.use = 123, verbose = FALSE)

# inspect stage, developmental age and dataset mixing
DimPlot(PVN_trajectory_dev, group.by = "trajectory_stage", pt.size = 0.1, shuffle = TRUE) +
  ggtitle("Trajectory stage") + 
  theme_void()

DimPlot(PVN_trajectory_dev, group.by = "gestational.week",
                       pt.size = 0.1, shuffle = TRUE) +
  ggtitle("Gestational week") + 
  theme_void()

DimPlot(PVN_trajectory_dev, group.by = "dataset",
                          pt.size = 0.1, shuffle = TRUE) +
  ggtitle("Dataset") + 
  theme_void()

# convert the Seurat object to a Monocle3 object
cds_pvn <- as.cell_data_set(PVN_trajectory_dev, assay = "RNA", reductions = "umap", default.reduction = "umap")

# identify Monocle clusters and partitions
set.seed(123)
cds_pvn <- cluster_cells(cds_pvn, reduction_method = "UMAP")

# learn separate trajectory graphs within each partition
cds_pvn <- learn_graph(cds_pvn, use_partition = TRUE)

# use the partition containing most GW6 Neuroepithelial/RG cells as root
early_cells_all <- which(colData(cds_pvn)$trajectory_stage == "Neuroepithelial/RG" & colData(cds_pvn)$gestational.week == "GW6")
main_partition <- names(which.max(table(partitions(cds_pvn)[early_cells_all])))
main_cells <- names(partitions(cds_pvn))[partitions(cds_pvn) == main_partition]
early_cells <- early_cells_all[partitions(cds_pvn)[early_cells_all] == main_partition]

# root the graph at the principal node occupied by most early cells
closest_vertex_all <- as.matrix(principal_graph_aux(cds_pvn)[["UMAP"]]$pr_graph_cell_proj_closest_vertex[colnames(cds_pvn), , drop = FALSE])
root_index <- as.integer(names(which.max(table(closest_vertex_all[early_cells, 1]))))
root_node <- igraph::V(principal_graph(cds_pvn)[["UMAP"]])$name[root_index]
cds_pvn <- order_cells(cds_pvn, root_pr_nodes = root_node)

# transfer pseudotime to Seurat and retain the main connected trajectory
pvn_pseudotime <- monocle3::pseudotime(cds_pvn)
PVN_trajectory_dev$pseudotime <- pvn_pseudotime[colnames(PVN_trajectory_dev)]
PVN_trajectory_dev$pseudotime[!is.finite(PVN_trajectory_dev$pseudotime)] <- NA
cds_pvn_main <- cds_pvn[, main_cells]
rowData(cds_pvn_main)$gene_short_name <- rownames(cds_pvn_main)

# extract the Monocle principal graph
pvn_graph <- principal_graph(cds_pvn_main)[["UMAP"]]

# assign each cell to its nearest principal-graph node
closest_vertex <- as.matrix(principal_graph_aux(cds_pvn_main)[["UMAP"]]$pr_graph_cell_proj_closest_vertex[colnames(cds_pvn_main), , drop = FALSE])
node_names <- V(pvn_graph)$name
node_index <- closest_vertex[, 1]

# convert node indices to graph-node names
closest_node <- if (all(as.character(node_index) %in% node_names)) {
  as.character(node_index)
  } else {
    node_names[as.integer(node_index)]
  }

names(closest_node) <- rownames(closest_vertex)
colData(cds_pvn_main)$closest_node <- closest_node[colnames(cds_pvn_main)]

cell_metadata <- as.data.frame(colData(cds_pvn_main)) %>%
  rownames_to_column("cell") %>%
  mutate(monocle_pseudotime = monocle3::pseudotime(cds_pvn_main)[cell])

# characterise four-edge neighbourhoods around all graph leaves
terminal_nodes <-V(pvn_graph)$name[igraph::degree(pvn_graph) == 1]
terminal_distances <- distances(pvn_graph, v = terminal_nodes, to = igraph::V(pvn_graph), weights = NA)

terminal_counts <- map_dfr(terminal_nodes, function(tip) {
  nearby_nodes <- colnames(terminal_distances)[terminal_distances[tip, ] <= 4]
  cell_metadata %>%
    filter(closest_node %in% nearby_nodes) %>%
    dplyr::count(trajectory_identity, name = "n") %>%
    mutate(terminal_node = tip, .before = 1)
})

# select the leaves most strongly occupied by each mature endpoint identity
select_tip <- function(identity, exclude = character()) {
  terminal_counts %>%
    filter(trajectory_identity == identity,
           !terminal_node %in% exclude) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    pull(terminal_node)
}

avp_oxt_tip <- select_tip("PVN: AVP-OXT")
avp_crh_tip <- select_tip("PVN: AVP-CRH")
non_nen_tip <- select_tip("PVN: non_NEN", c(avp_oxt_tip, avp_crh_tip))

late_terminal_nodes <- c("Lineage A: OXT/AVP-OXT" = avp_oxt_tip, "Lineage B: TRH/AVP-CRH" = avp_crh_tip, "Lineage C: non_NEN terminal" = non_nen_tip)

# trace each endpoint to the root and identify lineage-specific graph nodes
endpoint_paths <- lapply(late_terminal_nodes, function(tip) {
  path <- igraph::shortest_paths(pvn_graph, from = root_node, to = tip, output = "vpath")$vpath[[1]]
  igraph::V(pvn_graph)[path]$name
})

branch_specific_nodes <- imap(endpoint_paths, function(path_nodes, endpoint) {
  other_nodes <- unique(unlist(endpoint_paths[names(endpoint_paths) != endpoint]))
  setdiff(path_nodes, other_nodes)
})

# label lineage-specific arms, the common trunk and unmodelled side paths
cell_metadata <- cell_metadata %>%
  mutate(trajectory_branch = case_when(
    closest_node %in% branch_specific_nodes[[1]] ~ names(branch_specific_nodes)[1],
    closest_node %in% branch_specific_nodes[[2]] ~ names(branch_specific_nodes)[2],
    closest_node %in% branch_specific_nodes[[3]] ~ names(branch_specific_nodes)[3],
    closest_node %in% unique(unlist(endpoint_paths)) ~ "Shared trunk",
    TRUE ~ NA_character_))

colData(cds_pvn_main)$trajectory_branch <- cell_metadata$trajectory_branch[
  match(colnames(cds_pvn_main), cell_metadata$cell)]
PVN_trajectory_dev$trajectory_branch <- cell_metadata$trajectory_branch[
  match(colnames(PVN_trajectory_dev), cell_metadata$cell)]

# cell-identity palette
identity_colours <- c("Neuroepithelial/RG"   = "#5FA9A2",
                      "Cycling progenitor"   = "#C5A04A",
                      "IPC"                  = "#DD805E",
                      "Immature neuron"      = "#858B91",
                      
                      # OXT/AVP-OXT-related identities: purple
                      "PVN: OXT"             = "#A995C9",
                      "PVN: AVP-OXT"         = "#765996",
                      
                      # TRH/AVP-CRH-related identities: blue
                      "PVN: TRH"             = "#587F9E",
                      "PVN: AVP-CRH"         = "#78A8D1",
                      
                      # non-NEN lineage: teal
                      "PVN: non_NEN"         = "#D78FA9")


# pseudotime gradient
pseudotime_colours <- c("#30123B", "#7358A4", "#B85F78", "#DD805E", "#F3D35A")

# create function to plot uMAPs
plot_graph <- function(colour_by, title) {
  
  plot_cells(cds_pvn_main,
             color_cells_by = colour_by,
             label_cell_groups = FALSE,
             label_leaves = FALSE,
             label_branch_points = FALSE,
             label_roots = FALSE,
             trajectory_graph_color = "grey25",
             trajectory_graph_segment_size = 0.5,
             cell_size = 0.5) +
    ggtitle(title) +
    theme_void() +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5, margin = margin(b = 12)),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.key.size = unit(3, "mm"),
          legend.spacing.y = unit(0, "mm"),
          legend.box.spacing = unit(0, "pt"))
}


# plot UMAPs by cell identity, GW, pseudotime and trajectory branch
p_graph_identity <- plot_graph("trajectory_identity", "PVN identities") +
  scale_colour_manual(values = identity_colours, name = NULL) + 
  guides(colour = guide_legend(title = "PVN identities", override.aes = list(size = 1.5)))

# pseudotime UMAP
p_graph_pseudotime <- plot_graph("pseudotime", "Pseudotime") +
  scale_colour_gradientn(colours = pseudotime_colours,  name = "Pseudotime", na.value = "grey90")

# crop plots to the range occupied by actual cells
umap_coordinates <- reducedDims(cds_pvn_main)$UMAP
x_limits <- range(umap_coordinates[, 1], na.rm = TRUE)
y_limits <- range(umap_coordinates[, 2], na.rm = TRUE)
p_graph_identity <- p_graph_identity + coord_cartesian(xlim = x_limits, ylim = y_limits, expand = FALSE)
p_graph_pseudotime <- p_graph_pseudotime + coord_cartesian(xlim = x_limits, ylim = y_limits, expand = FALSE)

# plot together
(p_graph_identity | plot_spacer() | p_graph_pseudotime) +
  plot_layout(widths = c(1, 0.1, 1),
              guides = "collect")

ggsave("fig5.5a_pseudotime_umap.tiff", height = 3, width = 6.27, dpi = 600, compression = "lzw")

# violin plot of pseudotime by trajectory identity
p_stage_pseudotime <- VlnPlot(PVN_trajectory_dev, features = "pseudotime",
                              group.by = "trajectory_stage", pt.size = 0) +
  NoLegend() + RotatedAxis()


ggsave(file.path(figure_dir, "trajectory_umap_QC.pdf"),
       p_umap_overview, width = 14, height = 6)
ggsave(file.path(figure_dir, "trajectory_graph_overview.pdf"),
       p_trajectory_overview, width = 15, height = 6)
ggsave(file.path(figure_dir, "trajectory_branches.pdf"),
       p_graph_branch, width = 7, height = 6)
ggsave(file.path(figure_dir, "trajectory_stage_pseudotime.pdf"),
       p_stage_pseudotime, width = 7, height = 5)

# retain TFs associated with the PVN trajectory
trajectory_tfs <- trajectory_genes_filtered %>%
  filter(gene %in% human_tfs) %>%
  arrange(desc(morans_I))

# reusable pseudotime-expression plotting function
plot_pseudotime_genes <- function(genes, ncol = 4) {
  genes <- intersect(genes, rownames(PVN_trajectory_dev))
  
  expression_data <- FetchData(
    PVN_trajectory_dev,
    vars = c("pseudotime", "trajectory_stage", genes),
    cells = colnames(cds_pvn_main), layer = "data") %>%
    rownames_to_column("cell") %>%
    filter(is.finite(pseudotime)) %>%
    pivot_longer(all_of(genes), names_to = "gene", values_to = "expression") %>%
    mutate(gene = factor(gene, levels = genes))
  
  ggplot(expression_data, aes(pseudotime, expression)) +
    geom_smooth(method = "gam", formula = y ~ s(x, k = 6),
                colour = "#A995C9", se = FALSE, linewidth = 0.6) +
    facet_wrap(~gene, scales = "free_y", ncol = ncol) +
    labs(x = "Pseudotime", y = "Normalised expression") +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 9),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 8),
          axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3))
}

# plot global maturation, targeted TF and stress/endocrine programmes
maturation_genes <- c("SOX2", "TOP2A", "ASCL1", "DCX", "SNAP25", "SIM1")
plot_pseudotime_genes(maturation_genes, ncol = 3)

ggsave("fig5.5b_lineage_markers_smooth.tiff", height = 3, width = 6.27, dpi = 600, compression = "lzw")
