# load required packages
library(here)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)

# set working directory
setwd(here("results", "CELLECT"))

# import MAGMA results
magma_res <- read_csv("results_magma/prioritization.csv")

# import LDSC results
ldsc_res <- read_csv("results_ldsc/prioritization.csv")

# combine MAGMA and LDSC results and perfrom BH FDR correction
cellect_res <- bind_rows(magma_res %>% mutate(method = "MAGMA"),
                         ldsc_res %>% mutate(method = "LDSC")) %>%
  group_by(method) %>%
  mutate(padj = p.adjust(pvalue, method = "BH"),
         significance = case_when(padj < 0.001 ~ "***",
                                  padj < 0.01  ~ "**",
                                  padj < 0.05  ~ "*",
                                  TRUE ~ "ns")) %>%
  ungroup()

# inspect significant enrichments
cellect_res %>%
  filter(padj < 0.05) %>%
  arrange(method, padj) %>%
  select(method, gwas, annotation, beta, pvalue, padj)

# define cell population and disorder order
cellect_NEN_order <- c("AVP_CRH", "AVP_OXT", "OXT", "TRH", "GHRH", "PVN_non_NEN", "ARC_non_NEN")

cellect_GWAS_order <- c("ASD_Grove_2019", "MDD_PGC_2025", "ANX_Storm_2026", "PTSD_Nievergelt_2024",
                        "ADHD_Demontis_2023", "SCZ_Trubetskoy_2022", "BD_OConnell_2025")

# create MAGMA and LDSC enrichment and significance matrices
cellect_matrices <- lapply(c("MAGMA", "LDSC"), function(method_name) {
  
  data <- cellect_res %>%
    filter(method == method_name) %>%
    mutate(logp = -log10(pvalue))
  
  list(values = data %>%
      select(gwas, annotation, logp) %>%
      pivot_wider(names_from = annotation, values_from = logp) %>%
      column_to_rownames("gwas") %>%
      as.matrix(),
    
    significance = data %>%
      select(gwas, annotation, significance) %>%
      mutate(significance = replace(significance, significance == "ns", "")) %>%
      pivot_wider(names_from = annotation, values_from = significance) %>%
      column_to_rownames("gwas") %>%
      as.matrix())
  })

names(cellect_matrices) <- c("MAGMA", "LDSC")

# arrange rows and columns and simplify disorder labels
for (method_name in names(cellect_matrices)) {
  
  for (matrix_type in c("values", "significance")) {
    cellect_matrices[[method_name]][[matrix_type]] <- cellect_matrices[[method_name]][[matrix_type]][cellect_GWAS_order, cellect_NEN_order, drop = FALSE]
  }
  
  rownames(cellect_matrices[[method_name]]$values) <- c("ASD", "MDD", "ANX", "PTSD", "ADHD", "SCZ", "BD")
  }

# shorten population labels in both heatmap matrices
cellect_short_labels <- c("AVP-\nCRH", "AVP-\nOXT", "OXT", "TRH", "GHRH", "PVN\nnon-\nNEN", "ARC\nnon-\nNEN")

for (method_name in names(cellect_matrices)) {
  colnames(cellect_matrices[[method_name]]$values) <- cellect_short_labels
}

# create common colour scale across both methods
cellect_max <- max(cellect_matrices$MAGMA$values,  cellect_matrices$LDSC$values, na.rm = TRUE)
cellect_heatmap_colours <- colorRamp2(c(0, cellect_max / 2, cellect_max), c("#E9DFFC", "#7656AB"))


ht_opt$DIMNAME_PADDING <- unit(3, "mm")

# MAGMA heatmap
cellect_magma_heatmap <- Heatmap(cellect_matrices$MAGMA$values,
                                 name = "-log10(p)",
                                 col = cellect_heatmap_colours,
                                 cluster_rows = FALSE,
                                 cluster_columns = FALSE,
                                 row_names_gp = gpar(fontsize = 7),
                                 row_names_centered = TRUE,
                                 column_names_gp = gpar(fontsize = 7),
                                 column_names_rot = 0,
                                 column_names_centered = TRUE,
                                 column_names_max_height = unit(20, "mm"),
                                 column_title = "MAGMA",
                                 column_title_gp = gpar(fontsize = 9, fontface = "bold"),
                                 height = unit(60, "mm"),
                                 width = unit(60, "mm"),
                                 cell_fun = function(j, i, x, y, width, height, fill) {
                                   if (cellect_matrices$MAGMA$significance[i, j] != "") {
                                     grid.text(cellect_matrices$MAGMA$significance[i, j], x, y,
                                               gp = gpar(fontsize = 7, fontface = "bold"))
                                     }
                                   },
                                 heatmap_legend_param = list(title_gp = gpar(fontsize = 7),
                                                             labels_gp = gpar(fontsize = 7),
                                                             grid_width = unit(3, "mm"),
                                                             grid_height = unit(5, "mm")))

# LDSC heatmap
cellect_ldsc_heatmap <- Heatmap(cellect_matrices$LDSC$values,
                                name = "-log10(p)",
                                col = cellect_heatmap_colours,
                                cluster_rows = FALSE,
                                cluster_columns = FALSE,
                                row_names_gp = gpar(fontsize = 7),
                                row_names_centered = TRUE,
                                column_names_gp = gpar(fontsize = 7),
                                column_names_rot = 0,
                                column_names_centered = TRUE,
                                column_names_max_height = unit(20, "mm"),
                                column_title = "LDSC",
                                column_title_gp = gpar(fontsize = 9, fontface = "bold"),
                                height = unit(60, "mm"),
                                width = unit(60, "mm"),
                                cell_fun = function(j, i, x, y, width, height, fill) {
                                  if (cellect_matrices$LDSC$significance[i, j] != "") {
                                    grid.text(cellect_matrices$LDSC$significance[i, j], x, y,
                                              gp = gpar(fontsize = 7, fontface = "bold"))
                                    }
                                  },
                                heatmap_legend_param = list(title_gp = gpar(fontsize = 7),
                                                            labels_gp = gpar(fontsize = 7),
                                                            grid_width = unit(3, "mm"),
                                                            grid_height = unit(5, "mm")))

# display the heatmaps together
pdf("CELLECT_MAGMA_LDSC_heatmap.pdf", height = 2.9, width = 6, useDingbats = FALSE)
draw(cellect_magma_heatmap + cellect_ldsc_heatmap, ht_gap = unit(5, "mm"),  merge_legends = TRUE)
dev.off()