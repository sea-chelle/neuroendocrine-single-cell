# load packages
library(here)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)

# set working directory
setwd(here("plots"))

# combine TFs and simplify the labels used in the figure
tf_heatmap_genes <- bind_rows(tf_divergent_genes, tf_shared_genes) %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(display_program = recode(terminal_program, 
                                  "Strong OXT/AVP-OXT-enriched" = "OXT/AVP-OXT-enriched",
                                  "Strong TRH/AVP-CRH-enriched" = "TRH/AVP-CRH-enriched",
                                  "Strong non-NEN-enriched" = "non-NEN-enriched"),
         display_program = factor(display_program, levels = tf_program_levels)) %>%
  arrange(display_program)

# retain genes present in the expression matrix
tf_heatmap_genes <- tf_heatmap_genes %>%
  filter(gene %in% rownames(terminal_mean_expression), !is.na(display_program))

# extract expression in the same order as the annotation table
tf_expression <- terminal_mean_expression[tf_heatmap_genes$gene, , drop = FALSE]

# scale each gene across the three terminal groups
tf_expression_z <- t(scale(t(tf_expression)))

# remove genes without expression variance
keep_tf_genes <- apply(tf_expression_z, 1, function(x) all(is.finite(x)))
tf_expression_z <- tf_expression_z[keep_tf_genes, , drop = FALSE]

# restore exact agreement between annotation and expression rows
tf_heatmap_genes <- tf_heatmap_genes[match(rownames(tf_expression_z), tf_heatmap_genes$gene), , drop = FALSE]

# programme colours
tf_programme_colours <- c("OXT/AVP-OXT-enriched" = "#A995C9",
                          "TRH/AVP-CRH-enriched" = "#78A8D1",
                          "Shared NEN-enriched" = "#D78FA9",
                          "non-NEN-enriched" = "#5FA9A2",
                          "Shared across all terminals" = "#C5A95A")

# row annotation
tf_row_annotation <- rowAnnotation(Programme = tf_heatmap_genes$display_program,
                                   col = list(Programme = tf_programme_colours),
                                   show_annotation_name = FALSE,
                                   simple_anno_size = unit(3, "mm"),
                                   annotation_legend_param = list(title_gp = gpar(fontsize = 7),
                                                                  labels_gp = gpar(fontsize = 7),
                                                                  grid_width = unit(3, "mm"),
                                                                  grid_height = unit(3, "mm")))

# build TF heatmap
tf_heatmap <- Heatmap(tf_expression_z,
                      name = "Z-score",
                      col = colorRamp2(c(-1.2, 0, 1.2), c("#2674B8", "white", "#D9364A")),
                      left_annotation = tf_row_annotation,
                      row_split = tf_heatmap_genes$display_program,
                      cluster_rows = TRUE,
                      show_row_dend = FALSE,
                      cluster_row_slices = FALSE,
                      cluster_columns = FALSE,
                      row_title = NULL,
                      row_gap = unit(1.5, "mm"),
                      row_names_side = "right",
                      row_names_gp = gpar(fontsize = 7),
                      column_names_centered = TRUE,
                      column_names_side = "bottom",
                      column_names_rot = 0,
                      column_names_gp = gpar(fontsize = 8),
                      column_title = "Terminal PVN transcription factor programmes",
                      column_title_gp = gpar(fontsize = 8, fontface = "bold"),
                      heatmap_legend_param = list(title_gp = gpar(fontsize = 7),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_width = unit(3, "mm"),
                                                  grid_height = unit(3, "mm"),
                                                  at = c(-1.2, 0, 1.2)))

draw(tf_heatmap, merge_legends = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")

# programme labels displayed in the non-TF heatmap
non_tf_program_levels <- c("OXT/AVP-OXT-enriched", "TRH/AVP-CRH-enriched", "Shared NEN-enriched", "non-NEN-enriched",
                           "Shared across all terminals", "Defining neuropeptides")

# retain relevant programmes and create simplified display labels
non_tf_candidates <- non_tf_results %>%
  filter(terminal_program %in% c("Strong OXT/AVP-OXT-enriched", "Strong TRH/AVP-CRH-enriched", "Strong non-NEN-enriched",
                                 "Shared NEN-enriched", "non-NEN-enriched", "Shared across all terminals")) %>%
  mutate(display_program = recode(terminal_program,
                                  "Strong OXT/AVP-OXT-enriched" = "OXT/AVP-OXT-enriched",
                                  "Strong TRH/AVP-CRH-enriched" = "TRH/AVP-CRH-enriched",
                                  "Strong non-NEN-enriched" = "non-NEN-enriched"),
         plot_rank = case_when(display_program == "Shared across all terminals" ~ min_terminal_detection, TRUE ~ max_abs_logFC))

# select ten representative genes per displayed programme
top_non_tf_genes <- non_tf_candidates %>%
  group_by(display_program) %>%
  slice_max(plot_rank, n = 10, with_ties = FALSE) %>%
  ungroup()

# add defining neuropeptides explicitly
target_neuropeptides <- non_tf_results %>%
  filter(gene %in% c("OXT", "AVP", "CRH", "TRH")) %>%
  mutate(display_program = "Defining neuropeptides", plot_rank = Inf)

# combine unbiased and targeted genes
non_tf_heatmap_genes <- bind_rows(target_neuropeptides, top_non_tf_genes) %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(display_program = factor(display_program, levels = non_tf_program_levels)) %>%
  arrange(display_program) %>%
  filter(gene %in% rownames(terminal_mean_expression), !is.na(display_program))

# extract terminal mean expression
non_tf_expression <- terminal_mean_expression[non_tf_heatmap_genes$gene, , drop = FALSE]

# scale each gene across terminal groups
non_tf_expression_z <- t(scale(t(non_tf_expression)))

# remove genes with no expression variance
keep_non_tf_genes <- apply(non_tf_expression_z, 1, function(x) all(is.finite(x)))
non_tf_expression_z <- non_tf_expression_z[keep_non_tf_genes, , drop = FALSE]

# ensure annotation and expression matrices have identical row order
non_tf_heatmap_genes <- non_tf_heatmap_genes[match(rownames(non_tf_expression_z), non_tf_heatmap_genes$gene), , drop = FALSE]

# programme colours
non_tf_programme_colours <- c("OXT/AVP-OXT-enriched" = "#A995C9",
                              "TRH/AVP-CRH-enriched" = "#78A8D1",
                              "Shared NEN-enriched" = "#5FA9A2",
                              "non-NEN-enriched" = "#D78FA9",
                              "Shared across all terminals" = "#C5A95A",
                              "Defining neuropeptides" = "#DD805E")

# row annotation
non_tf_row_annotation <- rowAnnotation(Programme = non_tf_heatmap_genes$display_program,
                                       col = list(Programme = non_tf_programme_colours),
                                       show_annotation_name = FALSE,
                                       simple_anno_size = unit(3, "mm"),
                                       annotation_legend_param = list(title_gp = gpar(fontsize = 8),
                                                                      labels_gp = gpar(fontsize = 8),
                                                                      grid_width = unit(3, "mm"),
                                                                      grid_height = unit(3, "mm")))

# build heatmap
non_tf_heatmap <- Heatmap(non_tf_expression_z, 
                          name = "Z-score",
                          col = colorRamp2(c(-1, 0, 1), c("#2674B8", "white", "#D9364A")),
                          left_annotation = non_tf_row_annotation,
                          row_split = non_tf_heatmap_genes$display_program,
                          
                          # cluster genes within each programme without showing dendrograms
                          cluster_rows = TRUE,
                          show_row_dend = FALSE,
                          cluster_row_slices = FALSE,
                          cluster_columns = FALSE,
                          row_title = NULL,
                          row_gap = unit(1.5, "mm"),
                          row_names_side = "right",
                          row_names_gp = gpar(fontsize = 7),
                          column_names_side = "bottom",
                          column_names_rot = 0,
                          column_names_centered = TRUE,
                          column_names_gp = gpar(fontsize = 8),
                          column_title = "Terminal PVN protein-coding gene programmes",
                          column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                          heatmap_legend_param = list(title_gp = gpar(fontsize = 7),
                                                      labels_gp = gpar(fontsize = 7),
                                                      grid_width = unit(3, "mm"),
                                                      grid_height = unit(3, "mm"),
                                                      at = c(-1, 0, 1)))

draw(non_tf_heatmap, merge_legends = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")

# log-normalise the complete tradeSeq count matrix
trade_library_sizes <- colSums(counts_trade)
trade_logcounts <- counts_trade %*%
  Diagonal(x = 10000 / trade_library_sizes)
trade_logcounts@x <- log1p(trade_logcounts@x)
dimnames(trade_logcounts) <- dimnames(counts_trade)

# create function for plotting genes by psuedotime
plot_lineage_genes <- function(genes, plot_title = NULL, facet_ncol = 3) 
{ 
  # retain genes present in the expression matrix
  genes <- intersect(genes, rownames(trade_logcounts))
  plot_data <- map_dfr(seq_len(base::ncol(cell_weights)),
                       function(lineage_index) {
                         lineage_cells <- which(cell_weights[, lineage_index] > 0 & is.finite(normalised_pseudotime[, lineage_index]))
                         
                         # extract these before mutate to avoid name collisions
                         lineage_weights <- cell_weights[lineage_cells, lineage_index]
                         lineage_times <- normalised_pseudotime[lineage_cells, lineage_index]
                         expression_data <- as.data.frame(as.matrix(t(trade_logcounts[genes, lineage_cells, drop = FALSE])))
                         expression_data %>%
                           mutate(cell = rownames(cell_weights)[lineage_cells],
                                  normalised_pseudotime = lineage_times,
                                  lineage = lineage_names[lineage_index],
                                  lineage_weight = lineage_weights) %>%
                           pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expression")
                       })
  
  # preserve requested gene order
  plot_data$gene <- factor(plot_data$gene, levels = genes)
  
  # preserve legend order
  plot_data$lineage <- factor(plot_data$lineage, levels = lineage_names)
  
  # create plot
  ggplot(plot_data, aes(x = normalised_pseudotime, y = expression, colour = lineage)) +
    geom_smooth(aes(weight = lineage_weight), method = "gam", formula = y ~ s(x, k = 6), se = FALSE, linewidth = 0.6) +
    facet_wrap(~ gene, scales = "free_y", ncol = facet_ncol) +
    scale_colour_manual(values = c("OXT/AVP-OXT" = "#A995C9", "TRH/AVP-CRH" = "#78A8D1", "non-NEN" = "#D78FA9")) +
    labs(x = "Relative pseudotime", y = "Log-normalised expression", colour = "Lineage", title = plot_title) +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 9),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 9),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 8),
          legend.position = "bottom",
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.margin = margin(b = 0),
          legend.box.spacing = unit(0, "mm"),
          axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3))
}

# plot final selected TFs
headline_tf_genes <- c("NPAS4", "NR3C1", "NR3C2", "ISL1", "BCL6", "GTF2I")
headline_tf_plot <- plot_lineage_genes(genes = headline_tf_genes, facet_ncol = 3)
headline_tf_plot