# loading required packages
library(here)
library(Seurat)
library(tidyverse)
library(EWCE)
library(patchwork)

# set working directory to figures
setwd(here("plots"))

# subset NENs and corresponding PVN and ARC non-NEN populations
genetics_dev <- subset(nuclei_dev, subset = NEN_type == "NEN" |
                         (NEN_type != "NEN" & NEN_region %in% c("PVN", "ARC")))

# annotate NEN subtypes and regional non-NEN populations
genetics_dev$genetics_group <- ifelse(genetics_dev$NEN_type == "NEN", 
                                      as.character(genetics_dev$NEN_subtype),
                                      paste0(genetics_dev$NEN_region, " non-NEN"))

genetics_dev$genetics_group <- factor(genetics_dev$genetics_group,
                                      levels = c("AVP-CRH", "AVP-OXT", "OXT", "TRH", "PVN non-NEN", "GHRH", "ARC non-NEN"))

# check population sizes
table(genetics_dev$genetics_group, useNA = "ifany")

# import SFARI ASD-associated genes
SFARI <- read.csv(here("csv", "SFARI-Gene_genes_05-01-2026release_05-20-2026export.csv"))

# filter for category 1 non-syndromic genes
SFARI_cat1_genes <- SFARI %>%
  filter(gene.score == 1, syndromic == 0) %>%
  pull(gene.symbol) %>%
  unique()

# import MDD-associated genes from Open Targets and DisGeNET
OpenTarget <- read_tsv(here("csv", "opentarget_MDD.tsv"), show_col_types = FALSE)
DisGeNET <- read.delim(unz(here("csv", "DisGeNET_MDD.zip"), "search_result_C1269683.tsv"))
MDD_genes <- union(OpenTarget %>% filter(globalScore > 0.5) %>% pull(symbol), DisGeNET %>% filter(score > 0.5) %>% pull(Gene))

# calculate average expression of ASD-associated genes across neuronal populations
SFARI_avg_expr <- AverageExpression(genetics_dev,
                                    features = SFARI_cat1_genes,
                                    assays = "RNA",
                                    group.by = "genetics_group")$RNA

# scale expression across populations for each gene
SFARI_expr_scaled <- t(scale(t(as.matrix(SFARI_avg_expr))))

# remove genes with no variation in expression
SFARI_expr_scaled <- SFARI_expr_scaled[rowSums(is.finite(SFARI_expr_scaled)) == ncol(SFARI_expr_scaled), , drop = FALSE]


# calculate average expression of MDD-associated genes across neuronal populations
MDD_avg_expr <- AverageExpression(genetics_dev,
                                  features = MDD_genes,
                                  assays = "RNA",
                                  group.by = "genetics_group")$RNA

# scale expression across populations for each gene
MDD_expr_scaled <- t(scale(t(as.matrix(MDD_avg_expr))))

# remove genes with no variation in expression
MDD_expr_scaled <- MDD_expr_scaled[rowSums(is.finite(MDD_expr_scaled)) == ncol(MDD_expr_scaled), , drop = FALSE]

# plot MDD-associated gene expression
Heatmap(MDD_expr_scaled,
        name = "Z-score",
        col = colorRamp2(c(-2, 0, 2.5), c("#2674B8", "white", "#D9364A")),
        cluster_rows = TRUE,
        show_row_dend = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        column_names_gp = gpar(fontsize = 8),
        column_title = "MDD-associated genes")

# plot ASD-associated gene expression
Heatmap(SFARI_expr_scaled,
        name = "Z-score",
        col = colorRamp2(c(-2, 0, 2.5), c("#2674B8", "white", "#D9364A")),
        cluster_rows = TRUE,
        show_row_dend = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        column_names_gp = gpar(fontsize = 8),
        column_title = "SFARI category 1 ASD-associated genes")

# count and rank genes with highest average expression in each population
genetics_gene_counts <- bind_rows(tibble(population = colnames(SFARI_avg_expr)[max.col(as.matrix(SFARI_avg_expr), ties.method = "first")],
                                         disorder = "ASD"),
                                  tibble(population = colnames(MDD_avg_expr)[max.col(as.matrix(MDD_avg_expr), ties.method = "first")], 
                                         disorder = "MDD")) %>%
  count(disorder, population, name = "n_genes") %>%
  group_by(disorder) %>%
  arrange(desc(n_genes), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# print ranked counts
print(genetics_gene_counts, n = Inf)

# extract RNA expression and population annotations for EWCE
DefaultAssay(genetics_dev) <- "RNA"
EWCE_expr <- GetAssayData(genetics_dev, assay = "RNA", layer = "data")
EWCE_groups <- setNames(as.character(genetics_dev$genetics_group), colnames(genetics_dev))
EWCE_groups <- EWCE_groups[colnames(EWCE_expr)]
stopifnot(identical(colnames(EWCE_expr), names(EWCE_groups)), !anyNA(EWCE_groups))

# retain disease-associated genes represented in the expression data
EWCE_SFARI_genes <- intersect(SFARI_cat1_genes, rownames(EWCE_expr))
EWCE_MDD_genes <- intersect(MDD_genes, rownames(EWCE_expr))
cat("SFARI genes supplied/present:", length(SFARI_cat1_genes), "/", length(EWCE_SFARI_genes), "\n")
cat("MDD genes supplied/present:", length(MDD_genes), "/", length(EWCE_MDD_genes), "\n")

# generate EWCE reference for seven neuronal populations
EWCE_ctd_file <- generate_celltype_data(exp = EWCE_expr,
                                        annotLevels = list(genetics_group = EWCE_groups),
                                        groupName = "Developing_Hypothalamus",
                                        input_species = "human",
                                        output_species = "human",
                                        specificity_quantiles = TRUE,
                                        savePath = getwd())
EWCE_ctd <- load_rdata(EWCE_ctd_file)

# test ASD gene-set enrichment with gene-length and GC-content matching
set.seed(123)
EWCE_SFARI <- bootstrap_enrichment_test(sct_data = EWCE_ctd,
                                        hits = EWCE_SFARI_genes,
                                        bg = rownames(EWCE_ctd[[1]]$specificity),
                                        genelistSpecies = "human",
                                        sctSpecies = "human",
                                        annotLevel = 1,
                                        reps = 100000,
                                        geneSizeControl = TRUE,
                                        mtc_method = "BH")

# test MDD gene-set enrichment with the same settings
set.seed(123)
EWCE_MDD <- bootstrap_enrichment_test(sct_data = EWCE_ctd,
                                      hits = EWCE_MDD_genes,
                                      bg = rownames(EWCE_ctd[[1]]$specificity),
                                      genelistSpecies = "human",
                                      sctSpecies = "human",
                                      annotLevel = 1,
                                      reps = 100000,
                                      geneSizeControl = TRUE,
                                      mtc_method = "BH")

# inspect EWCE results
EWCE_SFARI$results
EWCE_MDD$results

# prepare EWCE fold-enrichment results for plotting
EWCE_plot <- bind_rows(EWCE_SFARI$results %>% rownames_to_column("population") %>% mutate(disorder = "ASD"),
                       EWCE_MDD$results %>% rownames_to_column("population") %>% mutate(disorder = "MDD")) %>%
  mutate(population = gsub("-non-NEN", " non-NEN", gsub("_", "-", population)), significant = q < 0.05)

# plot ASD EWCE fold enrichment
ASD_EWCE_plot <- EWCE_plot %>%
  filter(disorder == "ASD") %>%
  ggplot(aes(x = reorder(population, fold_change), y = fold_change, fill = significant)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  coord_flip(ylim = c(0, 1.5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("FALSE" = "grey80", "TRUE" = "#A995C9")) +
  labs(title = "ASD", x = NULL, y = "EWCE fold enrichment") +
  theme_classic() +
  theme(plot.title = element_text(size = 9, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"),
        legend.spacing.y = unit(0, "mm"),
        legend.box.spacing = unit(0, "pt"))

# plot MDD EWCE fold enrichment
MDD_EWCE_plot <- EWCE_plot %>%
  filter(disorder == "MDD") %>%
  ggplot(aes(x = reorder(population, fold_change), y = fold_change, fill = significant)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  coord_flip(ylim = c(0, 1.5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("FALSE" = "grey80", "TRUE" = "#A995C9")) +
  labs(title = "MDD", x = NULL, y = "EWCE fold enrichment") +
  theme_classic() +
  theme(plot.title = element_text(size = 9, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"),
        legend.spacing.y = unit(0, "mm"),
        legend.box.spacing = unit(0, "pt"))

# combine ASD and MDD EWCE plots with shared legend
(ASD_EWCE_plot | MDD_EWCE_plot) + 
  plot_layout(guides = "collect") + 
  theme(legend.position = "right")

# export combined figure
ggsave("EWCE_enrichment_combined.tiff", width = 6.5, height = 2, dpi = 600, compression = "lzw")

# calculate detection percentages for genes in ASD and MDD lists
genetics_pct_expr <- DotPlot(genetics_dev,
                             features = unique(c(EWCE_SFARI_genes, EWCE_MDD_genes)),
                             group.by = "genetics_group")$data %>%
  transmute(gene = features.plot, population = as.character(id), pct_expr = pct.exp)

# rank ASD-associated genes by specificity in selected populations
SFARI_specificity_ranked <- EWCE_ctd[[1]]$specificity[EWCE_SFARI_genes, c("AVP-CRH", "TRH", "GHRH"), drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "population", values_to = "specificity") %>%
  left_join(genetics_pct_expr, by = c("gene", "population")) %>%
  filter(pct_expr >= 10) %>%
  group_by(population) %>%
  arrange(desc(specificity), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# rank MDD-associated genes by specificity in selected populations
MDD_specificity_ranked <- EWCE_ctd[[1]]$specificity[EWCE_MDD_genes, c("AVP-CRH", "AVP-OXT", "OXT"), drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "population", values_to = "specificity") %>%
  left_join(genetics_pct_expr, by = c("gene", "population")) %>%
  filter(pct_expr >= 10) %>%
  group_by(population) %>%
  arrange(desc(specificity), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# retain top 20 detected genes per population for candidate selection
SFARI_top20 <- SFARI_specificity_ranked %>% filter(rank <= 20)
MDD_top20 <- MDD_specificity_ranked %>% filter(rank <= 20)
SFARI_top20 %>% print(n = Inf)
MDD_top20 %>% print(n = Inf)

# select ASD-associated candidate genes from the top 20 lists
SFARI_candidate_genes <- c("DSCAM", "NLGN4X", "NRXN1", "NRXN3", "SHANK2", "FOXP2", "AUTS2", "RUNX1T1", "GRIN2B")

# select MDD-associated candidate genes from the top 20 lists
MDD_candidate_genes <- c("NR3C1", "NR3C2", "CRH", "HTR2C", "NTRK2", "NEGR1", "MAOA", "GABRA2", "OPRM1")

# plot ASD-associated candidate genes
ASD_genes_dotplot <- DotPlot(genetics_dev,
                             features = SFARI_candidate_genes,
                             group.by = "genetics_group",
                             idents = c("AVP-CRH", "TRH", "GHRH"), 
                             dot.scale = 5) +
  scale_colour_gradient2(low = "#E3DAF0", high = "#7656AB", limits = c(-1.5, 1.5), breaks = c(-1.5, 0, 1.5)) +
  scale_size(limits = c(0, 100),  breaks = c(25, 50, 75), name = "Percent Expressed") +
  RotatedAxis() +
  labs(title = "ASD", x = NULL, y = NULL) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 9),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.position = "right",
        legend.key.size = unit(3, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.box.spacing = unit(0, "pt"),
        panel.spacing = unit(0.5, "lines"),
        panel.heights = unit(3 * 1, "cm"))

# plot MDD-associated candidate genes
MDD_genes_dotplot <- DotPlot(genetics_dev, features = MDD_candidate_genes, 
                             group.by = "genetics_group", 
                             idents = c("AVP-CRH", "AVP-OXT", "OXT"),
                             dot.scale = 5) +
  scale_colour_gradient2(low = "#E3DAF0", high = "#7656AB", limits = c(-1.5, 1.5), breaks = c(-1.5, 0, 1.5)) +
  scale_size(limits = c(0, 100),  breaks = c(25, 50, 75), name = "Percent Expressed") +
  RotatedAxis() +
  labs(title = "MDD", x = NULL, y = NULL) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 9),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.position = "right",
        legend.key.size = unit(3, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.box.spacing = unit(0, "pt"),
        panel.spacing = unit(0.5, "lines"),
        panel.heights = unit(3 * 1, "cm"))

# combine ASD and MDD EWCE plots with shared legend
(ASD_genes_dotplot + MDD_genes_dotplot) + 
  plot_layout(guides = "collect", widths = c(1, 1))

# export combined figure
ggsave("ASD_MDD_gene_expression.tiff", width = 6.5, height = 2.5, dpi = 600, compression = "lzw")
