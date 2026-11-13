# load required packages
library(here)
library(tidyverse)
library(EnhancedVolcano)

# set working directory 
setwd(here("plots"))

# select comparisons for volcano plots
volcano_results <- list(NEN_identity = NEN_DE_results$Overall, 
                        Sex_differences = NEN_DE_results$Male_vs_Female,
                        Development = dev_adult_DE_results$Overall)

# define genes to label on volcano plots
volcano_labels <- list(NEN_identity = c("SCGN", "OXT", "TH", "POU3F2", "AVP", "OTP", "PCSK1", "VGF", "GAL",  "ADRA1A", 
                                        "ADRA1B", "CADPS2", "GHRH", "BCL6", "CRH", "TENM2", "NXPH1", "AR", "HTR4", "DRD3", "SGK1"),
                       
                       Sex_differences = c("UTY", "NLGN4Y", "HTR2C", "ESR1","BAIAP3", "MAGEL2", "IGSF1",
                                           "XIST", "CALCB", "NR4A3", "GPR39", "DLL1"),
                       
                       Development = c("THRB", "FKBP5", "KCNQ3", "NEAT1", "OXR1", "PCSK1", "SNAP25", "SOX2", "ASCL1",
                                       "NKX2-1", "TUBB3", "GPC2", "PNOC", "IGFBP2"))

# define volcano plot titles
volcano_titles <- c(NEN_identity = "NEN vs non-NEN", 
                    Sex_differences = "Male vs female NEN", 
                    Development = "Adult vs developmental NEN")

# function to generate volcano plots
plot_NEN_volcano <- function(res, title, label_genes) {
  
  # create labels for selected genes
  labels <- ifelse(res$gene %in% label_genes, res$gene, "")
  
  # plot volcano plot
  EnhancedVolcano(res, x = "log2FoldChange", y = "padj",
                  lab = labels,
                  pCutoff = 0.05,
                  FCcutoff = 0.5,
                  cutoffLineWidth = 0.3,
                  labSize = 2,
                  pointSize = 1,
                  shape = 16,
                  colAlpha = 0.8,
                  legendIconSize = 2,
                  col = c("#b3b3b3", "#A995C9", "#DD805E", "#5FA9A2"),
                  subtitle = NULL,
                  drawConnectors = TRUE,
                  widthConnectors = 0.2,
                  colConnectors = "grey50",
                  arrowheads = FALSE,
                  max.overlaps = Inf) +
    ggtitle(title) +
    theme_bw() +
    theme(plot.caption = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 8),
          axis.title = element_text(size = 7),
          axis.text = element_text(size = 7),
          axis.line = element_line(linewidth = 0.3),
          legend.position = "none")
}

# generate volcano plots
NEN_volcano_plots <- lapply(names(volcano_results), function(name) {
  plot_NEN_volcano(res = volcano_results[[name]],
                   title = volcano_titles[[name]],
                   label_genes = volcano_labels[[name]])
})

names(NEN_volcano_plots) <- names(volcano_results)

# display volcano plots
NEN_volcano_plots$NEN_identity
NEN_volcano_plots$Sex_differences
NEN_volcano_plots$Development

# save volcano plots
for (name in names(NEN_volcano_plots)) {
  ggsave(plot = NEN_volcano_plots[[name]],
         filename = paste0(name, "_volcano.tiff"),
         height = 3.1, width = 3.1, dpi = 600, compression = "lzw")
}