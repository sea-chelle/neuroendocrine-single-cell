# load required library
library(here)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(igraph)
library(ggraph)
library(lme4)
library(lmerTest)

# set working directory 
setwd(here("plots"))

# specify directory for pySCENIC results
scenic_output_dir <- here("results", "pySCENIC", "output", "adult")

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

# top 20 regulons by RSS for each NEN subtype
scenic_mature_rss_balanced %>%
  dplyr::select(regulon, CRH, AVP_PVN, AVP_SON, OXT, TRH, GHRH) %>%
  pivot_longer(-regulon, names_to = "subtype", values_to = "RSS") %>%
  group_by(subtype) %>%
  slice_max(RSS, n = 20, with_ties = FALSE) %>%
  arrange(subtype, desc(RSS)) %>%
  print(n = Inf)

# select subtype-associated regulons to plot
scenic_subtype_selected <- list("CRH" = c("JUNB(+)", "ATF4(+)", "THRA(+)", "CEBPB(+)"),
                                "AVP_PVN" = c("CLOCK(+)", "NKX2-2(+)", "NR2F2(+)", "OTP(+)"),
                                "AVP_SON" = c("FOS(+)", "EGR1(+)", "MEIS1(+)", "CUX1(+)", "JUN(+)"),
                                "OXT" = c( "FOXP2(+)", "POU3F2(+)", "THAP11(+)", "NR2F1(+)"),
                                "TRH" = c("NFIX(+)", "NFIA(+)", "NHLH1(+)", "VAX1(+)", "DBP(+)"),
                                "GHRH" = c( "DLX1(+)", "FOXO1(+)", "NR3C1(+)", "GSX1(+)", 'BCL6(+)'))

# define subtype order for the mature-neuron heatmaps
scenic_subtype_order <- c("CRH", "AVP_PVN", "AVP_SON", "OXT", "TRH", "GHRH")

# extract balanced RSS values and convert to a regulon-by-subtype matrix
scenic_subtype_matrix <- scenic_mature_rss_balanced %>%
  dplyr::select(regulon, all_of(scenic_subtype_order)) %>%
  column_to_rownames("regulon") %>% as.matrix()
scenic_subtype_matrix <- scenic_subtype_matrix[unlist(scenic_subtype_selected), , drop = FALSE]

# calculate Z-score
scenic_subtype_scaled <- t(scale(t(scenic_subtype_matrix)))

# remove (+) from TF names
rownames(scenic_subtype_scaled) <- gsub("\\(\\+\\)", "", rownames(scenic_subtype_scaled))

# RSS heatmap
pdf("scenic_NEN_heatmap.pdf", height = 4, width = 3.5, useDingbats = FALSE)
Heatmap(scenic_subtype_scaled,
        name = "RSS",
        col = colorRamp2(c(-2, 0, 2), c("#2674B8", "white", "#D9364A")),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 7), 
        column_names_gp = gpar(fontsize = 7),
        column_names_rot = 0,
        column_names_centered = TRUE,
        column_labels = c("CRH", "PVN\nAVP", "SON\nAVP", "OXT", "TRH", "GHRH"),
        heatmap_legend_param = list(title_gp = gpar(fontsize = 7),
                                    labels_gp = gpar(fontsize = 7),
                                    grid_width = unit(3, "mm"),
                                    grid_height = unit(4, "mm"),
                                    at = c(-2, 0, 2)))
dev.off()

# choose regulators to illustrate subtype-specific activity
scenic_focus_regulons <- c("JUNB", "CLOCK", "FOS", "POU3F2", "NFIA", "NR3C1")

# reshape the unscaled, per-cell AUCell scores into long format for PVN NENs only
scenic_focus_AUC <- scenic_auc_annotated %>%
  filter(scenic_group %in% c("PVN: non-NEN", "ARC: non-NEN", "CRH", "AVP_PVN", "AVP_SON", "OXT", "TRH", "GHRH")) %>%
  dplyr::select(scenic_group, all_of(paste0(scenic_focus_regulons, "(+)"))) %>%
  pivot_longer(-scenic_group, names_to = "regulon", values_to = "AUC") %>%
  mutate(regulon = factor(sub("\\(\\+\\)$", "", regulon), levels = scenic_focus_regulons),
         scenic_group = factor(scenic_group,
                               levels = c("PVN: non-NEN", "ARC: non-NEN", "CRH", "AVP_PVN", "AVP_SON", "OXT", "TRH", "GHRH"),
                               labels = c("PVN non-NEN", "ARC non-NEN", "CRH", "PVN AVP", "SON AVP", "OXT", "TRH", "GHRH")))

# show the single-cell activity distributions as violins with inset boxplots
ggplot(scenic_focus_AUC, aes(x = scenic_group, y = AUC, fill = scenic_group)) +
  geom_violin(trim = TRUE, scale = "width", colour = "grey40", linewidth = 0.25) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.7, linewidth = 0.25) +
  facet_wrap(~regulon, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("SON AVP" = "#5FA9A2",
                               "PVN non-NEN" = "#B85F5A",
                               "GHRH" = "#D78FA9",
                               "OXT" = "#78A8D1",
                               "PVN AVP" = "#829B73",
                               "TRH" = "#A995C9",
                               "CRH" = "#C5A04A",
                               "ARC non-NEN" = "#DD805E")) +
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
ggsave("SCENIC_aucell_violin.tiff", height = 4.5, width = 5.5, dpi = 600, compression = "lzw")

scenic_regulon %>%
  filter(TF == "NR3C1") %>%
  dplyr::select(TF, TargetGenes)

# fix column names
scenic_regulon <- scenic_regulon %>%
  rename(TF = `...1`, TargetGenes = `Enrichment...9`)

# extract targeets for NR3C1
nr3c1_targets <- scenic_regulon %>%
  filter(TF == "NR3C1") %>%
  mutate(target = str_extract_all(TargetGenes, "(?<=\\(')[^']+(?=',)")) %>%
  dplyr::select(TF, target) %>%
  unnest(target) %>%
  distinct(TF, target)

nrow(nr3c1_targets)

# identify NR3C1 targets in the PVN TRH/AVP-CRH terminal programme
nr3c1_targets <- intersect(scenic_target_sets[["NR3C1"]], scenic_tradeSeq_sets[["TRH/AVP-CRH enriched"]])

# manually select biologically relevant targets for a readable network
# these are illustrative functional categories, not formal GO terms
nr3c1_genes <- list("Glucocorticoid" = c("FKBP5", "FKBP4", "KLF9", "NCOA2"),
                    "Neuronal excitability" = c("CACNB2", "GRIA2", "GABBR1"),
                    "Synaptic" = c("NLGN4X", "NLGN1", "NRXN3"),
                    "Endocrine" = c("IGSF1", "GNAS", "LEPR", "NPY1R", "MRAP2", "GHR", "ESR2", "THRB", "AGTR1"),
                    "TF" = c("FOXO1", "GSX1", "FOXP2", "VAX1", "NR2F2", "BCL6"))

# retain only selected genes verified within the NR3C1 overlap
nr3c1_network <- tibble(target = unlist(nr3c1_genes, use.names = FALSE),
                        category = rep(names(nr3c1_genes), lengths(nr3c1_genes))) %>% filter(target %in% nr3c1_targets$target)

# build the inferred NR3C1-to-target network; no target–target edges
nr3c1_graph <- graph_from_data_frame(nr3c1_network %>% transmute(from = "NR3C1", to = target),
                                     directed = FALSE,
                                     vertices = bind_rows(tibble(name = "NR3C1", category = "NR3C1"),
                                                          nr3c1_network %>% transmute(name = target, category)))

# place target genes in five functional arcs, keeping NR3C1 at the centre
nr3c1_layout <- tibble(name = V(nr3c1_graph)$name,
                       category = V(nr3c1_graph)$category) %>%
  filter(name != "NR3C1") %>%
  mutate(category = factor(category, levels = c("Glucocorticoid", "Neuronal excitability", "Synaptic", "Endocrine", "TF"))) %>%
  arrange(category) %>%
  mutate(angle = seq(0, 360 - 360 / n(), length.out = n())) %>%
  ungroup() %>%
  mutate(x = 2.5 * cos(angle * pi / 180), y = 2.5 * sin(angle * pi / 180)) %>%
  dplyr::select(name, category, x, y) %>%
  bind_rows(tibble(name = "NR3C1", category = "NR3C1", x = 0, y = 0)) %>%
  arrange(match(name, V(nr3c1_graph)$name))

# plot the network 
ggraph(nr3c1_graph, layout = "manual", x = nr3c1_layout$x, y = nr3c1_layout$y) +
  geom_edge_link(alpha = 0.5, colour = "grey60", linewidth = 0.3) +
  geom_node_point(aes(colour = category, size = name == "NR3C1")) +
  geom_node_text(aes(label = ifelse(name == "NR3C1", "", name)), repel = TRUE, size = 2, box.padding = 0.2, point.padding = 0.6,
                 fontface = "italic", max.overlaps = Inf) +
  geom_node_text(aes(label = ifelse(name == "NR3C1", "NR3C1", "")),
                 nudge_y = 0.5, size = 2, fontface = "bold.italic") +
  scale_colour_manual(values = c("NR3C1" = "#A995C9",
                                 "Glucocorticoid" = "#D78FA9",
                                 "Neuronal excitability" = "#78A8D1",
                                 "Synaptic" = "#C5A04A",
                                 "Endocrine" = "#5FA9A2",
                                 "TF" = "#DD805E"),
                      breaks = c("Glucocorticoid", "Neuronal excitability", "Synaptic", "Endocrine", "TF"),
                      name = "Functional category") +
  guides(colour = guide_legend(nrow = 3, byrow = TRUE)) +
  scale_size_manual(values = c("FALSE" = 2.5, "TRUE" = 7), guide = "none") +
  coord_equal(clip = "off") + 
  theme_void(base_size = 8) +
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.spacing.y = unit(0, "mm"),
        legend.key.height = unit(3, "mm"),
        legend.key.spacing.y = unit(0, "mm"))

ggsave("SCENIC_NR3C1_network_adult.tiff", height = 3.3, width = 3.3, dpi = 600, compression = "lzw")

# identify all regulon columns
scenic_regulon_cols <- grep("\\(\\+\\)$", colnames(scenic_auc_sex), value = TRUE)

# calculate donor-level mean AUCell activity for all regulons
scenic_donor_AUC <- scenic_auc_sex %>%
  dplyr::select(cell, all_of(scenic_regulon_cols)) %>%
  left_join(scenic_metadata %>%
              dplyr::select(cell, donor, sex, NEN_status), by = "cell") %>%
  group_by(donor, sex, NEN_status) %>%
  summarise(across(all_of(scenic_regulon_cols), mean, na.rm = TRUE), n_cells = n(), .groups = "drop") %>%
  pivot_longer(cols = all_of(scenic_regulon_cols), names_to = "regulon", values_to = "AUC") %>%
  mutate(sex = factor(sex, levels = c("Female", "Male")), NEN_status = factor(NEN_status, levels = c("non_NEN", "NEN")))

n_distinct(scenic_donor_AUC$regulon)

# fit sex × NEN models across all regulons
scenic_sex_models <- scenic_donor_AUC %>%
  group_by(regulon) %>%
  group_modify(~ {
    model <- tryCatch(lmer(AUC ~ sex * NEN_status + (1 | donor), data = .x), error = function(e) NULL)
    
    if (is.null(model)) {
      return(tibble(effect = c("sex", "NEN_status", "sex:NEN_status"), pvalue = NA_real_, singular = NA))
    }
    
    as.data.frame(anova(model, type = 3)) %>%
      tibble::rownames_to_column("effect") %>%
      transmute(effect, pvalue = `Pr(>F)`, singular = isSingular(model))
  }) %>%
  ungroup()

# adjust p-values separately for each model effect
scenic_sex_models <- scenic_sex_models %>%
  group_by(effect) %>%
  mutate(FDR = p.adjust(pvalue, method = "BH")) %>%
  ungroup()

# count significant regulons
scenic_sex_models %>%
  group_by(effect) %>%
  summarise(n_tested = sum(!is.na(pvalue)), 
            n_FDR05 = sum(FDR < 0.05, na.rm = TRUE), 
            n_singular = sum(singular, na.rm = TRUE),
            .groups = "drop")

# sex main effect
scenic_sex_models %>%
  filter(effect == "sex") %>%
  arrange(FDR) %>%
  print(n = 50)

# NEN main effect
scenic_sex_models %>%
  filter(effect == "NEN_status") %>%
  arrange(FDR) %>%
  print(n = 20)

# sex × NEN interaction
scenic_sex_models %>%
  filter(effect == "sex:NEN_status") %>%
  arrange(FDR) %>%
  print(n = 50)

# define NEN and non-NEN populations
scenic_metadata_sex <- scenic_metadata %>%
  mutate(NEN_type = if_else(scenic_group %in% c("AVP_PVN", "AVP_SON", "CRH", "OXT", "TRH", "GHRH"), "NEN", "non-NEN"),
         sex_group = paste(sex, NEN_type, sep = "_"))

# format data for plotting
scenic_auc_sex <- scenic_auc %>%
  left_join(scenic_metadata_sex %>%
              dplyr::select(cell, sex, NEN_type, sex_group), by = "cell")

# choose regulators for sex- and NEN-associated activity
scenic_sex_regulons <- c("NR3C1", "FOS", "MEIS1", "STAT3", "PGR", "EGR1")

# reshape unscaled per-cell AUCell scores
scenic_sex_AUC <- scenic_auc_sex %>%
  dplyr::select(sex_group, all_of(paste0(scenic_sex_regulons, "(+)"))) %>%
  pivot_longer(-sex_group, names_to = "regulon", values_to = "AUC") %>%
  mutate(regulon = factor(sub("\\(\\+\\)$", "", regulon), levels = scenic_sex_regulons),
         sex_group = factor(sex_group,
                            levels = c("Male_NEN", "Male_non-NEN", "Female_NEN", "Female_non-NEN"),
                            labels = c("Male NEN", "Male non-NEN", "Female NEN", "Female non-NEN")))

# show single-cell regulon activity distributions
ggplot(scenic_sex_AUC, aes(x = sex_group, y = AUC, fill = sex_group)) +
  geom_violin(trim = TRUE, scale = "width", colour = "grey40", linewidth = 0.25) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white",
               alpha = 0.7, linewidth = 0.25) +
  facet_wrap(~regulon, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Male NEN" = "#78A8D1",
                               "Male non-NEN" = "#47749B",
                               "Female non-NEN" = "#B65679",
                               "Female NEN" = "#D78FA9")) +
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

ggsave("SCENIC_sex_aucell_violin_adult.tiff", height = 4, width = 4, dpi = 600, compression = "lzw")
