# load required packages
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

# define eligble protein-coding genes
go_results_table <- tradeSeq_results %>%
  filter(gene_class %in% c("TF", "non-TF"), !technical_gene, !str_detect(gene, fixed("///")))

# define dynamic developmental programmes
dynamic_go_sets <- list("Shared developmental dynamics" = go_results_table %>%
                          filter(is.finite(association_FDR), is.finite(pattern_FDR),
                                 association_FDR < 0.001, pattern_FDR >= 0.05) %>%
                          pull(gene) %>% unique(),
                        
                        "Lineage-divergent dynamics" = go_results_table %>%
                          filter(is.finite(pattern_FDR), pattern_FDR < 0.001, pattern_fcMedian >= 1) %>%
                          pull(gene) %>% unique())

# define terminal programmes using symmetrical criteria
terminal_go_sets <- list("OXT/AVP-OXT enriched" = go_results_table %>% 
                           filter(is.finite(endpoint_FDR),
                                  endpoint_FDR < 0.001,
                                  detection_OXT_AVP_OXT >= 5,
                                  logFC1_2 >= 1,
                                  logFC1_3 >= 1) %>%
                           pull(gene) %>%
                           unique(),
                         
                         "TRH/AVP-CRH enriched" = go_results_table %>%
                           filter(is.finite(endpoint_FDR),
                                  endpoint_FDR < 0.001,
                                  detection_TRH_AVP_CRH >= 5,
                                  logFC1_2 <= -1,
                                  logFC2_3 >= 1) %>%
                           pull(gene) %>%
                           unique(),
                         
                         "Shared NEN enriched" = go_results_table %>%
                           filter(is.finite(endpoint_FDR), 
                                  endpoint_FDR < 0.001,
                                  
                                  # expressed in both NEN terminal populations
                                  detection_OXT_AVP_OXT >= 5,
                                  detection_TRH_AVP_CRH >= 5,
                                  
                                  # similar between the two NEN branches
                                  abs(logFC1_2) < 1,
                                  
                                  # enriched in both NEN branches relative to non-NEN
                                  logFC1_3 >= 1,
                                  logFC2_3 >= 1) %>%
                           pull(gene) %>%
                           unique(),
                         
                         "non-NEN enriched" = go_results_table %>%
                           filter(is.finite(endpoint_FDR),
                                  endpoint_FDR < 0.001,
                                  detection_non_NEN >= 5,
                                  logFC1_3 <= -1,
                                  logFC2_3 <= -1) %>%
                           pull(gene) %>%
                           unique())

# check gene set sizes
bind_rows(tibble(analysis = "Developmental dynamics",
                 programme = names(dynamic_go_sets),
                 n_genes = lengths(dynamic_go_sets)),
           tibble(analysis = "Terminal identity",
                  programme = names(terminal_go_sets),
                  n_genes = lengths(terminal_go_sets)))

# background for terminal-programme GO
terminal_go_background <- go_results_table %>%
  filter(is.finite(endpoint_pvalue)) %>%
  pull(gene) %>%
  unique()

# GO Biological Process enrichment for terminal programmes
terminal_go_results <- compareCluster(geneCluster = terminal_go_sets,
                                      fun = "enrichGO",
                                      OrgDb = org.Hs.eg.db,
                                      keyType = "SYMBOL",
                                      ont = "BP",
                                      universe = terminal_go_background,
                                      pAdjustMethod = "BH",
                                      pvalueCutoff = 0.05,
                                      qvalueCutoff = 0.05,
                                      minGSSize = 10,
                                      maxGSSize = 500)

# convert GO results as dataframe
terminal_go_df <- as.data.frame(terminal_go_results)
terminal_go_df %>%
  group_by(Cluster) %>%
  slice_min(p.adjust, n = 10, with_ties = FALSE) %>%
  select(Cluster, ID, Description, GeneRatio, p.adjust)

# print all significant terminal GO terms
for (programme in unique(terminal_go_df$Cluster)) {
  
  cat("\n\n", programme, "\n", strrep("-", nchar(as.character(programme))), "\n", sep = "")
  terminal_go_df %>%
    filter(Cluster == programme, p.adjust < 0.05) %>%
    arrange(p.adjust) %>%
    pull(Description) %>%
    print(quote = FALSE)
  }

# select terms to plot
terminal_selected_terms <- tribble(~Cluster,                   ~Description,                                      ~PlotLabel,
                                   # OXT/AVP-OXT
                                   "OXT/AVP-OXT enriched",    "mitochondrial translation",                       "Mitochondrial translation",
                                   "OXT/AVP-OXT enriched",    "oxidative phosphorylation",                       "Oxidative phosphorylation",
                                   "OXT/AVP-OXT enriched",    "ribonucleotide biosynthetic process",             "Ribonucleotide biosynthesis",
                                   "OXT/AVP-OXT enriched",    "protein folding",                                 "Protein folding",
                                   
                                   # TRH/AVP-CRH
                                   "TRH/AVP-CRH enriched",    "modulation of chemical synaptic transmission",     "Synaptic transmission",
                                   "TRH/AVP-CRH enriched",    "regulation of membrane potential",                 "Membrane potential",
                                   "TRH/AVP-CRH enriched",    "neuropeptide signaling pathway",                   "Neuropeptide signalling",
                                   "TRH/AVP-CRH enriched",    "circadian behavior",                               "Circadian behaviour",
  
                                   # shared NEN
                                   "Shared NEN enriched",     "response to steroid hormone",                      "Steroid hormone response",
                                   "Shared NEN enriched",     "regulation of synaptic plasticity",                "Synaptic plasticity",
                                   "Shared NEN enriched",     "neurotransmitter uptake",                          "Neurotransmitter uptake",
                                   "Shared NEN enriched",     "long-term synaptic potentiation",                  "Long-term potentiation",
                                   
                                   # non-NEN
                                   "non-NEN enriched",        "asymmetric cell division",                         "Asymmetric cell division",
                                   "non-NEN enriched",        "neural precursor cell proliferation",              "Neural precursor proliferation",
                                   "non-NEN enriched",        "central nervous system neuron differentiation",    "CNS neuron differentiation",
                                   "non-NEN enriched",        "forebrain development",                            "Forebrain development")

# extract the selected enrichment results
terminal_go_plot_data <- terminal_go_df %>%
  inner_join(terminal_selected_terms, by = c("Cluster", "Description")) %>%
  mutate(GeneRatio_numeric = vapply(strsplit(GeneRatio, "/"),
                                    function(x) {
                                      as.numeric(x[1]) / as.numeric(x[2])
                                      },
                                    numeric(1)),
         neglog10FDR = -log10(p.adjust),
         Cluster = recode(as.character(Cluster), 
                          "OXT/AVP-OXT enriched" = "OXT/AVP-OXT",
                          "TRH/AVP-CRH enriched" = "TRH/AVP-CRH",
                          "Shared NEN enriched" = "NEN",
                          "non-NEN enriched" = "non-NEN"),
         Cluster = factor(Cluster, levels = c("OXT/AVP-OXT", "TRH/AVP-CRH", "NEN", "non-NEN")),
         PlotLabel = factor(PlotLabel, levels = rev(unique(terminal_selected_terms$PlotLabel))))

# plot terminal GO programmes
ggplot(terminal_go_plot_data, aes(x = GeneRatio_numeric, y = PlotLabel, size = Count, colour = neglog10FDR)) +
  geom_point() +
  facet_wrap(~Cluster, scales = "free_y", ncol = 2) +
  scale_colour_gradient(low = "#EEEAF7", high = "#7656AB") +
  scale_size_continuous(range = c(1, 5)) +
  scale_x_continuous(breaks = seq(0.02, 0.08, 0.02), limits = c(0.01, 0.09), expand = expansion(add = c(0.005, 0.005))) +
  labs(x = "Gene ratio", y = NULL, size = "Count", colour = expression(-log[10]("FDR"))) +
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
        panel.spacing = unit(0.5, "lines"),
        panel.heights = unit(10 * 0.5, "cm"))
ggsave("fig5.6c_tradeSeq_GO.tiff", height = 2.5, width = 6.27, dpi = 600, compression = "lzw")