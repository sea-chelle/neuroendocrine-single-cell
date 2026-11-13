# load required packages
library(here)
library(tidyverse)
library(Seurat)

# set working directory
setwd(here("plots"))

# subset all adult neurons
neuron_adult <- subset(seurat_adult, subset = NEN_type != "non_neuron")
Idents(neuron_adult) <- "NEN_subtype"

# identify markers for each NEN subtype against all other neurons
neuron_markers <- FindAllMarkers(neuron_adult,
                                 assay = "RNA",
                                 only.pos = TRUE,
                                 min.pct = 0.25,
                                 logfc.threshold = 0.25)

# filter for strongly enriched and widely detected neuronal markers
neuron_markers_filtered <- neuron_markers %>%
  mutate(pct_diff = pct.1 - pct.2) %>%
  filter(p_val_adj < 0.05, avg_log2FC >= 1, pct.1 >= 0.7, pct_diff >= 0.3)

# calculate average expression for filtered markers
neuron_avg_expression <- AverageExpression(neuron_adult,
                                           assays = "RNA",
                                           features = unique(neuron_markers_filtered$gene),
                                           group.by = "NEN_subtype",
                                           layer = "data")$RNA

# add average expression and rank markers
neuron_markers_filtered <- neuron_markers_filtered %>%
  mutate(avg_expression = mapply(function(gene, cluster) neuron_avg_expression[gene, gsub("_", "-", as.character(cluster))],
                                 gene, cluster)) %>%
  arrange(cluster, desc(avg_log2FC))

# filter for sufficiently expressed genes
neuron_markers_filtered <- neuron_markers_filtered %>%
  filter(avg_expression >= 1) %>%
  arrange(cluster, desc(avg_log2FC))

# view candidate merker genes  
neuron_markers_filtered %>%
    as_tibble() %>%
    print(n = Inf, width = Inf)

# define adult NEN co-markers
adult_NEN_comarkers <- list(CRH = c("CARTPT", "RASD1", "PENK", "RAMP1", "NR3C1"),
                            AVP_PVN = c("THSD4", "RMST"),
                            AVP_SON = c("PCSK1", "ADRA1A", "ADRA1B", "CADPS2", "KCNMB4"),
                            OXT = c("SCGN", "PTGFR", "TH", "CADPS2", "ADRA1A"),
                            TRH = c("BMPR1B", "NFIA", "NFIB", "TRHDE", "HTR2C", "GHR"),
                            GHRH = c("GAL", "DLK1", "ECEL1", "FOXP2", "GHR"))

# combine markers without duplicates and reverse order for plotting
adult_NEN_marker_plot <- rev(unique(unlist(adult_NEN_comarkers)))

# set order for plotting
NEN_adult$NEN_subtype <- factor(NEN_adult$NEN_subtype, levels = adult_NEN_order)

# plot dotplot of putative markers
DotPlot(NEN_adult, assay = "RNA", group.by = "NEN_subtype",
        features = adult_NEN_marker_plot,
        dot.scale = 3,
        cols = c("grey100", "#7656A8")) +
  coord_flip() +
  guides(colour = guide_colourbar(title = "Average\nexpression"), size = guide_legend(title = "Percent\nexpressed")) +
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
        panel.heights = unit(25 * 0.3, "cm"),
        panel.widths = unit(length(adult_NEN_order) * 0.5, "cm"))

ggsave("NEN_marker_dotplot.tiff", height = 3.5, width = 2.3, dpi = 600, compression = "lzw")

# subset all AVP neurons
AVP_adult <- subset(NEN_adult, subset = NEN_subtype %in% c("AVP_PVN", "AVP_SON"))
Idents(AVP_adult) <- "NEN_subtype"

# identify markers for each AVP PVN vs SON neurons
AVP_markers <- FindMarkers(AVP_adult,
                           ident.1 = "AVP_PVN",
                           ident.2 = "AVP_SON",
                           assay = "RNA",
                           min.pct = 0.25,
                           logfc.threshold = 0.25)

# add gene names and population labels
AVP_markers <- AVP_markers %>%
  rownames_to_column("gene") %>%
  mutate(cluster = ifelse(avg_log2FC > 0, "AVP_PVN", "AVP_SON")) %>%
  arrange(p_val_adj)

# calculate average expression for AVP markers
AVP_avg_expression <- AverageExpression(AVP_adult,
                                        assays = "RNA",
                                        features = AVP_markers$gene,
                                        group.by = "NEN_subtype",
                                        layer = "data")$RNA

# add average expression for the enriched AVP population
AVP_markers <- AVP_markers %>%
  mutate(avg_expression = AVP_avg_expression[cbind(match(gene, rownames(AVP_avg_expression)),
                                                   match(gsub("_", "-", cluster), colnames(AVP_avg_expression)))])

# filter for strongly enriched markers for each AVP population
AVP_markers_filtered <- AVP_markers %>%
  mutate(pct_target = ifelse(cluster == "AVP_PVN", pct.1, pct.2),
         pct_other = ifelse(cluster == "AVP_PVN", pct.2, pct.1),
         pct_diff = pct_target - pct_other) %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) >= 1, pct_target >= 0.7, pct_diff >= 0.3, avg_expression >= 1) %>%
  arrange(cluster, desc(abs(avg_log2FC)))

# view candidate merker genes  
AVP_markers_filtered %>%
  as_tibble() %>%
  print(n = Inf, width = Inf)

# define AVP population-specific co-markers
AVP_comarkers <- list(AVP_PVN = c("NNAT", "SOX5", "TAFA2", "BAIAP3", "PCDH15"),
                      AVP_SON = c("ADRA1A", "RGS6", "PCSK1", "CADPS2", "KCNMB4"))

# combine markers and reverse order for plotting
AVP_marker_plot <- rev(unique(unlist(AVP_comarkers)))

# set plotting order
AVP_adult$NEN_subtype <- factor(AVP_adult$NEN_subtype, levels = c("AVP_PVN", "AVP_SON"))

# plot dotplot of putative AVP co-markers
DotPlot(AVP_adult, assay = "RNA", group.by = "NEN_subtype",
        features = AVP_marker_plot,
        dot.scale = 3,
        cols = c("grey100", "#7656A8")) +
  coord_flip() +
  guides(colour = guide_colourbar(title = "Average\nexpression"), size = guide_legend(title = "Percent\nexpressed")) +
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
        panel.heights = unit(10 * 0.4, "cm"),
        panel.widths = unit(length(adult_NEN_order) * 0.4, "cm"))

ggsave("AVP_marker_dotplot.tiff", height = 2.4, width = 2.1, dpi = 600, compression = "lzw")