# load required library
library(here)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(igraph)
library(ggraph)

# set working directory 
setwd(here("plots"))

# specify directory for pySCENIC results
scenic_output_dir <- here("results", "pySCENIC", "output")

# import per-cell AUCell regulon activity scores
scenic_auc <- read_csv(file.path(scenic_output_dir, "pyscenic_auc_matrix.csv.gz"), show_col_types = FALSE)

# import balanced RSS scores for NEN versus non-NEN comparisons
scenic_NEN_rss_balanced <- read_csv(file.path(scenic_output_dir, "pyscenic_NEN_RSS_balanced.csv"), show_col_types = FALSE)

# import balanced RSS scores for mature PVN and ARC subtypes
scenic_mature_rss_balanced <- read_csv(file.path(scenic_output_dir, "pyscenic_mature_RSS_balanced.csv"), show_col_types = FALSE)

# import pySCENIC motif enrichment results containing predicted TF–target pairs
scenic_regulon <- read_csv(file.path(scenic_output_dir, "pyscenic_regulons.csv"), show_col_types = FALSE)

# add cell population annotations to per-cell AUCell scores
scenic_auc_annotated <- scenic_auc %>%
  left_join(scenic_metadata %>% 
              dplyr::select(cell, scenic_group), by = "cell")

# developmental regulon activity - combine the four PVN NEN subtypes into one developmental population and retain progenitors, immature neurons and other populations as separate groups
scenic_auc_development <- scenic_auc_annotated %>%
  mutate(developmental_group = case_when(scenic_group %in% c("PVN: OXT", "PVN: AVP-OXT", "PVN: AVP-CRH", "PVN: TRH") ~ "PVN NEN",
                                         scenic_group == "ARC: GHRH" ~ "ARC NEN",
                                         scenic_group == "PVN: non-NEN" ~ "PVN non-NEN",
                                         scenic_group == "ARC: non-NEN" ~ "ARC non-NEN",
                                         TRUE ~ scenic_group))

# calculate mean AUCell activity for each regulon and developmental population
scenic_development_master <- scenic_auc_development %>%
  select(-cell, -scenic_group) %>%
  group_by(developmental_group) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(-developmental_group, names_to = "regulon", values_to = "mean_auc") %>%
  pivot_wider(names_from = developmental_group, values_from = mean_auc)

# select representative TF regulons for the developmental heatmap
scenic_development_selected <- c("EMX2(+)", "SOX2(+)", "PAX6(+)", "HMGA1(+)", "FOXM1(+)", "HES1(+)", "ASCL1(+)", "TCF4(+)", "RFX3(+)", "ATF2(+)", "ATF4(+)", "CREB3(+)",
                                 "OTP(+)", "NEUROD1(+)", "EBF1(+)", "MESP1(+)", "NR2F6(+)", "THRA(+)", "KLF5(+)", "THRB(+)", "PBX1(+)", "GSX1(+)", "HMX2(+)", "DLX6(+)")

# define the order of cell populations along development
scenic_development_order <- c("Neuroepithelial/RG", "Cycling progenitor", "IPC", "Immature neuron", "PVN NEN", "PVN non-NEN", "ARC NEN", "ARC non-NEN")

# convert regulon activity into a numeric matrix
scenic_development_matrix_selected <- scenic_development_master %>%
  column_to_rownames("regulon") %>% as.matrix()

scenic_development_matrix_selected <- scenic_development_matrix_selected[scenic_development_selected, scenic_development_order, drop = FALSE]

# Z-score mean AUCell activity separately for each regulon (row)
scenic_development_scaled_selected <- t(scale(t(scenic_development_matrix_selected)))

# plot the developmental regulon activity heatmap
pdf("scenic_dev_heatmap.pdf", height = 5, width = 6.27, useDingbats = FALSE)
Heatmap(scenic_development_scaled_selected,
        name = "Mean AUCell\nZ-score",
        col = colorRamp2(c(-1.5, 0, 2.5), c("#2674B8", "white", "#D9364A")),
        cluster_rows = FALSE, cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8), 
        column_names_rot = 90,
        heatmap_legend_param = list(title_gp = gpar(fontsize = 7),
                                    labels_gp = gpar(fontsize = 7),
                                    grid_width = unit(3, "mm"),
                                    grid_height = unit(4, "mm"),
                                    at = c(-1.5, 0, 2.5)))
dev.off()

# shared NEN versus non-NEN regulon specificity (PVN/ARC)
# classify regulons by whether RSS favours NEN or non-NEN populations
# classification requires the same direction of specificity in PVN and ARC
scenic_NEN_master <- scenic_NEN_rss_balanced %>%
  mutate(identity_category = case_when(`PVN NEN` > `PVN non-NEN` & `ARC NEN` > `ARC non-NEN` ~ "Shared NEN",
                                       `PVN non-NEN` > `PVN NEN` & `ARC non-NEN` > `ARC NEN` ~ "Shared non-NEN",
                                       TRUE ~ "Other"))

# select shared and region-associated NEN/non-NEN regulons
scenic_NEN_selected <- list("Shared NEN" = c("DLX2(+)", "DLX5(+)", "EGR1(+)", "FOSB(+)"),
                            "Shared non-NEN" = c("EBF1(+)", "NFIA(+)", "PLAGL1(+)"),
                            "PVN NEN" = c("ONECUT2(+)", "NR2F6(+)", "JUNB(+)", "THRA(+)"),
                            "ARC NEN" = c("GSX1(+)", "HMX2(+)", "FOXO1(+)", "THRB(+)"),
                            "PVN non-NEN" = c("NR2F1(+)", "POU3F2(+)", "OTP(+)"),
                            "ARC non-NEN" = c("AR(+)", "NR3C1(+)", "CLOCK(+)", "FOXO3(+)"))

# prepare the balanced RSS matrix for selected neuronal populations
scenic_NEN_matrix <- scenic_NEN_master %>%
  select(regulon, `PVN NEN`, `ARC NEN`, `PVN non-NEN`, `ARC non-NEN`) %>%
  column_to_rownames("regulon") %>% as.matrix()
scenic_NEN_matrix <- scenic_NEN_matrix[unlist(scenic_NEN_selected), , drop = FALSE]

# Z-score RSS across populations for each regulon
scenic_NEN_scaled <- t(scale(t(scenic_NEN_matrix)))

# plot NEN versus non-NEN regulon specificity heatmap
pdf("scenic_NEN_heatmap.pdf", height = 5, width = 6.27, useDingbats = FALSE)
Heatmap(scenic_NEN_scaled,
        name = "RSS\nZ-score",
        col = colorRamp2(c(-2, 0, 2), c("#2674B8", "white", "#D9364A")),
        row_split = factor(rep(names(scenic_NEN_selected), lengths(scenic_NEN_selected)),
                     levels = names(scenic_NEN_selected)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 8), 
        column_names_gp = gpar(fontsize = 8),
        column_names_rot = 45,
        row_title_gp = gpar(fontsize = 8),
        row_title_rot = 0,
        heatmap_legend_param = list(title_gp = gpar(fontsize = 7),
                                    labels_gp = gpar(fontsize = 7),
                                    grid_width = unit(3, "mm"),
                                    grid_height = unit(4, "mm"),
                                    at = c(-2, 0, 2)))
dev.off()

# top 15 regulons by RSS for each NEN subtype
scenic_mature_rss_balanced %>%
  select(regulon, `PVN: AVP-CRH`, `PVN: AVP-OXT`,
         `PVN: OXT`, `PVN: TRH`, `ARC: GHRH`) %>%
  pivot_longer(-regulon, names_to = "subtype", values_to = "RSS") %>%
  group_by(subtype) %>%
  slice_max(RSS, n = 15, with_ties = FALSE) %>%
  arrange(subtype, desc(RSS)) %>%
  print(n = 75)

# select subtype-associated regulons to plot
scenic_subtype_selected <- list("AVP-CRH" = c("POU6F1(+)", "FOSL2(+)", "GMEB2(+)", "NR3C1(+)"),
                                "AVP-OXT" = c("IRF1(+)", "EGR1(+)", "JUNB(+)", "ATF3(+)", "THRA(+)"),
                                "OXT" = c("NR2F6(+)", "NR2F2(+)", "KMT2B(+)"),
                                "TRH" = c("UNCX(+)", "NR2F1(+)", "NEUROD1(+)", "OTP(+)"),
                                "GHRH" = c("GSX1(+)", "HMX2(+)", "FOXO1(+)", "DLX6(+)", "THRB(+)"))

# define subtype order for the mature-neuron heatmaps
scenic_subtype_order <- c("PVN: AVP-CRH", "PVN: AVP-OXT", "PVN: OXT", "PVN: TRH", "ARC: GHRH")

# extract balanced RSS values and convert to a regulon-by-subtype matrix
scenic_subtype_matrix <- scenic_mature_rss_balanced %>%
  select(regulon, all_of(scenic_subtype_order)) %>%
  column_to_rownames("regulon") %>% as.matrix()
scenic_subtype_matrix <- scenic_subtype_matrix[unlist(scenic_subtype_selected), , drop = FALSE]

# shorten subtype names for heatmap
colnames(scenic_subtype_matrix) <- sub("^PVN: |^ARC: ", "", colnames(scenic_subtype_matrix))
colnames(scenic_subtype_matrix) <- c("\nAVP-\nCRH", "\nAVP-\nOXT", "OXT", "TRH", "GHRH")

# calculate Z-score
scenic_subtype_scaled <- t(scale(t(scenic_subtype_matrix)))

# RSS heatmap
pdf("scenic_NEN_heatmap.pdf", height = 4, width = 3, useDingbats = FALSE)
Heatmap(scenic_subtype_scaled,
        name = "RSS",
        col = colorRamp2(c(-2, 0, 2), c("#2674B8", "white", "#D9364A")),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 7), 
        column_names_gp = gpar(fontsize = 7),
        column_names_rot = 0,
        column_names_centered = TRUE,
        heatmap_legend_param = list(title_gp = gpar(fontsize = 7),
                                    labels_gp = gpar(fontsize = 7),
                                    grid_width = unit(3, "mm"),
                                    grid_height = unit(4, "mm"),
                                    at = c(-2, 0, 2)))
dev.off()

# calculate mean per-cell AUCell activity for the selected regulons
scenic_subtype_activity <- scenic_auc_annotated %>%
  filter(scenic_group %in% scenic_subtype_order) %>%
  select(scenic_group, all_of(unlist(scenic_subtype_selected))) %>%
  group_by(scenic_group) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(-scenic_group, names_to = "regulon", values_to = "AUCell") %>%
  pivot_wider(names_from = scenic_group, values_from = AUCell)

# convert mean AUCell scores to a matrix in the desired regulon order
scenic_subtype_AUC_matrix <- scenic_subtype_activity %>%
  column_to_rownames("regulon") %>% as.matrix()
scenic_subtype_AUC_matrix <- scenic_subtype_AUC_matrix[unlist(scenic_subtype_selected), scenic_subtype_order, drop = FALSE]

# shorten subtype names for heatmap
colnames(scenic_subtype_AUC_matrix) <- sub("^PVN: |^ARC: ", "", colnames(scenic_subtype_AUC_matrix))
colnames(scenic_subtype_AUC_matrix) <- c("\nAVP-\nCRH", "\nAVP-\nOXT", "OXT", "TRH", "GHRH")

# Z-score AUCell activity across cell populations within each regulon
scenic_subtype_AUC_scaled <- t(scale(t(scenic_subtype_AUC_matrix)))

# plot mean AUCell activity for the selected mature neuronal subtypes
pdf("scenic_NEN_AUCell_heatmap.pdf", height = 5, width = 3, useDingbats = FALSE)
Heatmap(scenic_subtype_AUC_scaled,
        name = "Mean AUCell\nZ-score", 
        col = colorRamp2(c(-2, 0, 2), c("#2674B8", "white", "#D9364A")),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 7), 
        column_names_gp = gpar(fontsize = 7),
        column_names_rot = 0,
        column_names_centered = TRUE,
        heatmap_legend_param = list(title_gp = gpar(fontsize = 7),
                                    labels_gp = gpar(fontsize = 7),
                                    grid_width = unit(3, "mm"),
                                    grid_height = unit(4, "mm"),
                                    at = c(-2, 0, 2)))
dev.off()

# replace the two-line header with readable column names
colnames(scenic_regulon) <- c("TF", "MotifID", "AUC", "NES", "MotifSimilarityQvalue", "OrthologousIdentity", "Annotation", "Context", "TargetGenes", "RankAtMax")

# remove the two header rows from the imported CSV
scenic_regulon <- scenic_regulon %>% slice(-(1:2))

# extract gene symbols from Python-style TargetGenes lists and remove duplicate TF–target pairs
scenic_targets <- scenic_regulon %>%
  select(TF, TargetGenes) %>%
  filter(!is.na(TF), !is.na(TargetGenes)) %>%
  mutate(target = stringr::str_extract_all(TargetGenes, "(?<=\\(')[^']+(?=',)")) %>%
  select(TF, target) %>%
  unnest(target) %>%
  distinct(TF, target)

# retain eligible genes and remove technical annotations / ambiguous gene IDs
scenic_tradeSeq_genes <- tradeSeq_results %>%
  filter(gene_class %in% c("TF", "non-TF"),
         !coalesce(technical_gene, TRUE),
         !str_detect(gene, fixed("///")))

# define PVN developmental gene sets from the existing tradeSeq tests
scenic_tradeSeq_sets <- list("Shared developmental dynamics" = scenic_tradeSeq_genes %>%
                               filter(is.finite(association_FDR), is.finite(pattern_FDR),
                                      association_FDR < 0.001, pattern_FDR >= 0.05) %>% pull(gene) %>% unique(),
                             "Lineage-divergent dynamics" = scenic_tradeSeq_genes %>%
                               filter(is.finite(pattern_FDR), pattern_FDR < 0.001,
                                      pattern_fcMedian >= 1) %>% pull(gene) %>% unique(),
                             "OXT/AVP-OXT enriched" = scenic_tradeSeq_genes %>%
                               filter(is.finite(endpoint_FDR), endpoint_FDR < 0.001,
                                      detection_OXT_AVP_OXT >= 5, logFC1_2 >= 1, logFC1_3 >= 1) %>% pull(gene) %>% unique(),
                             "TRH/AVP-CRH enriched" = scenic_tradeSeq_genes %>%
                               filter(is.finite(endpoint_FDR), endpoint_FDR < 0.001,
                                      detection_TRH_AVP_CRH >= 5, logFC1_2 <= -1, logFC2_3 >= 1) %>% pull(gene) %>% unique(),
                             "Shared NEN enriched" = scenic_tradeSeq_genes %>%
                               filter(is.finite(endpoint_FDR), endpoint_FDR < 0.001,
                                      detection_OXT_AVP_OXT >= 5, detection_TRH_AVP_CRH >= 5,
                                      abs(logFC1_2) < 1, logFC1_3 >= 1, logFC2_3 >= 1) %>% pull(gene) %>% unique(),
                             "non-NEN enriched" = scenic_tradeSeq_genes %>% filter(is.finite(endpoint_FDR), endpoint_FDR < 0.001,
                                                                                   detection_non_NEN >= 5, logFC1_3 <= -1, logFC2_3 <= -1) %>% pull(gene) %>% unique())


# define the common gene universe used for the Fisher tests
scenic_tradeSeq_background <- intersect(scenic_tradeSeq_genes$gene, scenic_targets$target)

# restrict every developmental programme to that gene universe
scenic_tradeSeq_sets <- lapply(scenic_tradeSeq_sets, intersect, y = scenic_tradeSeq_background)

# fisher's exact test: all TFs x six PVN programmes
# collect the unique predicted targets of each TF, within the same background
scenic_target_sets <- split(scenic_targets$target, scenic_targets$TF)
scenic_target_sets <- lapply(scenic_target_sets, intersect, y = scenic_tradeSeq_background)

# test all TF × developmental-programme combinations for overrepresentation
scenic_tradeSeq_enrichment <- expand_grid(TF = names(scenic_target_sets),
                                          programme = names(scenic_tradeSeq_sets)) %>%
  rowwise() %>%
  mutate(n_targets = length(scenic_target_sets[[TF]]),
         n_programme = length(scenic_tradeSeq_sets[[programme]]),
         overlap = length(intersect(scenic_target_sets[[TF]], scenic_tradeSeq_sets[[programme]])),
   
          # expected overlap if genes were randomly associated within the background
         expected = n_targets * n_programme / length(scenic_tradeSeq_background),
         pvalue = fisher.test(matrix(c(overlap, n_targets - overlap,  n_programme - overlap, length(scenic_tradeSeq_background) - n_targets - n_programme + overlap), nrow = 2), alternative = "greater")$p.value) %>% 
  ungroup() %>%
  mutate(FDR = p.adjust(pvalue, method = "BH"),
         
    # enrichment strength = observed overlap / expected overlap
    if_else(expected > 0, overlap / expected, NA_real_),
    overlap_pct = if_else(n_targets > 0, 100 * overlap / n_targets, NA_real_)) %>%
  select(TF, programme, n_targets, n_programme, overlap,
         overlap_pct, expected, fold_enrichment, pvalue, FDR)

# compact summary of results table
scenic_tradeSeq_enrichment %>%
  filter(FDR < 0.05) %>% count(programme, name = "n_significant_regulons") %>%
  print(n = Inf)

# select TFs to plot
scenic_integration_selected <- list("Shared developmental dynamics" = "ZNF454",
                                    "Shared NEN" = c("JUNB", "THRA", "FOS", "EGR1", "JUN", "FOSB", "DLX2", "DLX5", "CEBPD", "TCF7L2", "ONECUT2"),
                                    "OXT/AVP-OXT" = c("MESP1", "ZNF580", "MZF1", "IRF1", "HMX2", "ATF3", "ILF2"),
                                    "TRH/AVP-CRH" = c("NR3C1", "THRB", "BCL6", "CLOCK", "FOSL2", "NFKB1"),
                                    "non-NEN" = c("NFIA", "PLAGL1", "POU3F2", "NR2F1", "OTP", "UNCX"))

# combine the list into unique TF names for plotting
scenic_plot_TFs <- unique(unlist(scenic_integration_selected, use.names = FALSE))

# set the biological ordering of the six developmental programmes
scenic_programme_order <- c("Shared developmental dynamics", "Lineage-divergent dynamics", "Shared NEN enriched", "OXT/AVP-OXT enriched", "TRH/AVP-CRH enriched", "non-NEN enriched")

# include only FDR-significant associations for the selected TFs
scenic_enrichment_plot <- scenic_tradeSeq_enrichment %>%
  filter(FDR < 0.05, TF %in% scenic_plot_TFs) %>%
  mutate(programme = factor(programme, levels = scenic_programme_order))

# plot
ggplot(scenic_enrichment_plot, aes(x = fold_enrichment, y = reorder(TF, fold_enrichment))) +
  geom_point(aes(size = overlap), colour = "#A995C9") +
  facet_wrap(~programme, scales = "free_y", ncol = 2,
    labeller = as_labeller(c("Shared developmental dynamics" = "Shared\ndynamics",
                             "Lineage-divergent dynamics" = "Divergent\ndynamics",
                             "Shared NEN enriched" = "Shared NEN",
                             "OXT/AVP-OXT enriched" = "OXT/AVP-OXT",
                             "TRH/AVP-CRH enriched" = "TRH/AVP-CRH",
                             "non-NEN enriched" = "non-NEN" ))) +
  scale_size_continuous(name = "Overlapping targets", range = c(1, 3)) +
  labs(x = "Fold enrichment", y = NULL) +
  theme_bw() +
  theme(panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.4),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(3, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.box.spacing = unit(0, "pt"),
        legend.position = "bottom")
ggsave("SCENIC_tradeSeq_dotplot.tiff", height = 6, width = 3.3, dpi = 600, compression = "lzw")

# choose regulators to illustrate subtype-specific activity
scenic_focus_regulons <- c("NR3C1", "FOSB", "EGR1", "FOSL2", "THRA", "THRB")

# reshape the unscaled, per-cell AUCell scores into long format for PVN NENs only
scenic_focus_AUC <- scenic_auc_annotated %>%
  filter(scenic_group %in% c("PVN: non-NEN", "PVN: OXT", "PVN: AVP-OXT", "PVN: AVP-CRH", "PVN: TRH")) %>%
  select(scenic_group, all_of(paste0(scenic_focus_regulons, "(+)"))) %>%
  pivot_longer(-scenic_group, names_to = "regulon", values_to = "AUC") %>%
  mutate(
    regulon = factor(sub("\\(\\+\\)$", "", regulon), levels = scenic_focus_regulons),
    scenic_group = factor(
      scenic_group,
      levels = c("PVN: non-NEN", "PVN: OXT", "PVN: AVP-OXT", "PVN: AVP-CRH", "PVN: TRH"),
      labels = c("non-NEN", "OXT", "AVP-OXT", "AVP-CRH", "TRH")
    )
  )

# show the single-cell activity distributions as violins with inset boxplots
ggplot(scenic_focus_AUC, aes(x = scenic_group, y = AUC, fill = scenic_group)) +
  geom_violin(trim = TRUE, scale = "width", colour = "grey40", linewidth = 0.25) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.7, linewidth = 0.25) +
  facet_wrap(~regulon, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("non-NEN" = "#5FA9A2",
                               "OXT" = "#D78FA9",
                               "AVP-OXT" = "#78A8D1",
                               "AVP-CRH" = "#A995C9",
                               "TRH" = "#C5A04A")) +
  labs(x = NULL, y = "AUCell activity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 7),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.title = element_text(size = 7),
        strip.text = element_text(face = "bold", size = 8),
        panel.spacing = unit(0.2, "lines"),
        legend.position = "none")
ggsave("SCENIC_aucell_violin.tiff", height = 4, width = 3.1, dpi = 600, compression = "lzw")

# identify NR3C1 targets in the PVN TRH/AVP-CRH terminal programme
nr3c1_targets <- intersect(scenic_target_sets[["NR3C1"]], scenic_tradeSeq_sets[["TRH/AVP-CRH enriched"]])

# manually select biologically relevant targets for a readable network
# these are illustrative functional categories, not formal GO terms
nr3c1_genes <- list("HPA" = c("CRH", "NR3C2", "RASD1", "HTR2C", "GPR158", "PDE4D"),
                    "Neuronal excitability" = c("CACNA1D", "CACNB2", "HCN1", "KCNMA1", "NALCN", "GABBR2", "GRM5"),
                    "Neurodevelopment" = c("NLGN4X", "NLGN1", "NRXN3", "DSCAM", "AHI1", "ROBO2", "SLIT3"),
                    "Endocrine" = c("IGSF1", "GNAS", "TACR3"))

# retain only selected genes verified within the NR3C1 overlap
nr3c1_network <- tibble(target = unlist(nr3c1_genes, use.names = FALSE),
                        category = rep(names(nr3c1_genes), lengths(nr3c1_genes))) %>% filter(target %in% nr3c1_targets)

# build the inferred NR3C1-to-target network; no target–target edges
nr3c1_graph <- graph_from_data_frame(nr3c1_network %>% transmute(from = "NR3C1", to = target),
                                     directed = FALSE,
                                     vertices = bind_rows(tibble(name = "NR3C1", category = "NR3C1"),
                                                          nr3c1_network %>% transmute(name = target, category)))

# place target genes in four functional arcs, keeping NR3C1 at the centre
nr3c1_layout <- tibble(name = V(nr3c1_graph)$name, category = V(nr3c1_graph)$category) %>%
  filter(name != "NR3C1") %>%
  mutate(start = case_when(category == "HPA" ~ 10,
                           category == "Neuronal excitability" ~ 100,
                           category == "Neurodevelopment" ~ 190,
                           category == "Endocrine" ~ 280)) %>%
  group_by(category) %>%
  mutate(angle = seq(first(start), first(start) + 70, length.out = n())) %>%
  ungroup() %>%
  mutate(x = 2.5 * cos(angle * pi / 180), y = 2.5 * sin(angle * pi / 180)) %>%
  select(name, category, x, y) %>%
  bind_rows(tibble(name = "NR3C1", category = "NR3C1", x = 0, y = 0)) %>%
  arrange(match(name, V(nr3c1_graph)$name))

# plot the network 
ggraph(nr3c1_graph, layout = "manual", x = nr3c1_layout$x, y = nr3c1_layout$y) +
  geom_edge_link(alpha = 0.5, colour = "grey60", linewidth = 0.3) +
  geom_node_point(aes(colour = category, size = name == "NR3C1")) +
  geom_node_text(aes(label = ifelse(name == "NR3C1", "", name)), repel = TRUE, size = 2, fontface = "italic", max.overlaps = Inf) +
  geom_node_text(aes(label = ifelse(name == "NR3C1", "NR3C1", "")),
                 nudge_y = 0.5, size = 2, fontface = "bold.italic") +
  scale_colour_manual(values = c("NR3C1" = "#A995C9",
                                 "HPA" = "#D78FA9",
                                 "Neuronal excitability" = "#78A8D1",
                                 "Neurodevelopment" = "#C5A04A",
                                 "Endocrine" = "#5FA9A2"),
                      breaks = c("HPA", "Neuronal excitability", "Neurodevelopment", "Endocrine"),
                      name = "Functional category") +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_size_manual(values = c("FALSE" = 2.5, "TRUE" = 7), guide = "none") +
  coord_equal(clip = "off") + 
  theme_void(base_size = 8) +
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.spacing.y = unit(0, "mm"),
        legend.key.height = unit(3, "mm"),
        legend.key.spacing.y = unit(0, "mm"))
ggsave("SCENIC_NR3C1_network.tiff", height = 2.2, width = 3, dpi = 600, compression = "lzw")

