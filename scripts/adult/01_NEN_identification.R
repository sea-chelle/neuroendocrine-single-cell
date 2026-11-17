# load required packages
library(here)
library(readxl)
library(tidyverse)
library(Seurat)

# set working directory
setwd(here("data", "tadross_2025"))

# import tadross 2025 seurat object (inc. 3 samples from herb/siletti dataset)
seurat_adult <- readRDS("human_HYPOMAP_snRNASeq.rds")

# load supplementary tables
excel_sheets("41586_2024_8504_MOESM3_ESM.xlsx")

# import cluster annotations and spatial assignments
adult_annotations <- read_excel("41586_2024_8504_MOESM3_ESM.xlsx", sheet = "6_cluster_tree_annotations")
adult_spatial <- read_excel("41586_2024_8504_MOESM3_ESM.xlsx", sheet = "12_region_abundance_assignment")

# inspect columns
colnames(adult_annotations)
colnames(adult_spatial)

# screen neuropeptide expression across neuronal C4 clusters
adult_NEN_screen <- DotPlot(subset(seurat_adult, subset = C0_named == "C0-3 Neuron"),
                            features = c("CRH", "AVP", "OXT", "TRH", "GHRH", "GNRH1", "SST", "TH", "DDC"),
                            group.by = "C4_named", assay = "RNA")$data

# define neuropeptides and anatomical regions of NENs
NEN_regions <- list(CRH = c("PVN"),
                    AVP = c("PVN", "SON"),
                    OXT = c("PVN", "SON"),
                    TRH = c("PVN"),
                    GHRH = c("ARC"),
                    GNRH1 = c("POA", "MPOA", "LPOA"),
                    SST = c("Periventricular"),
                    TH = c("ARC"),
                    DDC = c("ARC"))

# screen clusters with >=20% marker expression and compatible anatomical assignments
adult_NEN_candidates <- adult_NEN_screen %>%
  mutate(cluster_id = sub(" .*", "", as.character(id))) %>%
  left_join(adult_spatial %>%
              filter(grepl("^C4-", `C4 cluster`)) %>%
              select(`C4 cluster`, region_threshold_manual,
             `Region with highest abdunance score`),
            by = c("cluster_id" = "C4 cluster")) %>%
  filter(features.plot %in% names(NEN_regions), 
         pct.exp >= 20,
         mapply(function(marker, region) {
           is.na(region) || region == "NA" || region %in% NEN_regions[[as.character(marker)]]
      },
      features.plot, region_threshold_manual)) %>%
  arrange(features.plot, desc(avg.exp))

# number of candidates per neuropeptide
adult_NEN_candidates %>%
  count(features.plot)

# screen PVN/SON-associated neuropeptides
adult_PVN_candidates <- adult_NEN_candidates %>%
  filter(features.plot %in% c("CRH", "AVP", "OXT", "TRH"),
         `Region with highest abdunance score` %in% c("PVN", "SON"))

# rank neuropeptide expression
adult_PVN_candidates %>%
  select(cluster_id, features.plot, region_threshold_manual,
         `Region with highest abdunance score`, pct.exp, avg.exp) %>%
  arrange(features.plot, desc(avg.exp))

# define PVN neuroendocrine and preautonomic markers
PVN_NEN_markers <- c("CRH", "AVP", "OXT", "TRH", "AGTR1", "NR3C1", "SCGN", 
                     "CALB1", "KCNMB4", "FOXP1", "TAC1", "TH", "DDC", "NFIB", "THRB")

PVN_preautonomic_markers <- c("ADARB2", "RELN", "SNCA", "NTNG1", "CACNA1G")

# calculate neuroendocrine and preautonomic marker expression across candidate clusters
adult_PVN_expression <- DotPlot(subset(seurat_adult, 
                                       cells = colnames(seurat_adult)[seurat_adult$C4_named %in% unique(adult_PVN_candidates$id)]),
                                features = unique(c(PVN_NEN_markers, PVN_preautonomic_markers)),
                                group.by = "C4_named",
                                assay = "RNA")$data

# summarise co-marker expression for each cluster
adult_PVN_expression <- adult_PVN_expression %>%
  mutate(cluster_id = sub(" .*", "", as.character(id))) %>%
  select(cluster_id, marker = features.plot, pct.exp, avg.exp) %>%
  pivot_wider(names_from = marker, values_from = c(pct.exp, avg.exp),  names_glue = "{marker}_{.value}")

# examine CRH candidates and co-markers
adult_PVN_expression %>%
  filter(cluster_id %in% adult_PVN_candidates$cluster_id[adult_PVN_candidates$features.plot == "CRH"]) %>%
  select(cluster_id, starts_with("CRH_"),
         starts_with("AGTR1_"), starts_with("NR3C1_"), starts_with("SCGN_")) %>%
  arrange(desc(CRH_avg.exp)) 

# examine AVP candidates and co-markers
adult_PVN_expression %>%
  filter(cluster_id %in% adult_PVN_candidates$cluster_id[adult_PVN_candidates$features.plot == "AVP"]) %>%
  select(cluster_id, starts_with("AVP_"),
         starts_with("TH_"), starts_with("DDC"), starts_with("TAC1_"), starts_with("CALB1_"), starts_with("KCNMB4_"), starts_with("FOXP1_")) %>%
  arrange(desc(AVP_avg.exp)) %>%
  print(n = Inf, width = Inf)

# examine OXT candidates and co-markers
adult_PVN_expression %>%
  filter(cluster_id %in% adult_PVN_candidates$cluster_id[adult_PVN_candidates$features.plot == "OXT"]) %>%
  select(cluster_id, starts_with("OXT_"), starts_with("CALB1_"), starts_with("KCNMB4_"), starts_with("FOXP1_")) %>%
  arrange(desc(OXT_avg.exp))

# examine TRH candidates and co-markers
adult_PVN_expression %>%
  filter(cluster_id %in% adult_PVN_candidates$cluster_id[adult_PVN_candidates$features.plot == "TRH"]) %>%
  select(cluster_id, starts_with("TRH_"), starts_with("AGTR1_"), starts_with("NFIB_"), starts_with("THRB_")) %>%
  arrange(desc(TRH_avg.exp))

# select PVN/SON candidate populations following neuropeptide, co-marker and anatomical assessment
PVN_NEN_ids <- c("C4-293", "C4-295", "C4-296", "C4-297", "C4-298", "C4-312", "C4-313", "C4-314", "C4-315")

# subset shortlisted populations
PVN_NEN_subset <- subset(seurat_adult, cells = colnames(seurat_adult)[sub(" .*", "", seurat_adult$C4_named) %in% PVN_NEN_ids])
Idents(PVN_NEN_subset) <- "C4_named"

# examine neuroendocrine marker expression using stacked violins
VlnPlot(PVN_NEN_subset, features = PVN_NEN_markers, pt.size = 0, stack = TRUE, flip = TRUE)

# examine preautonomic marker expression using stacked violins
VlnPlot(PVN_NEN_subset, features = PVN_preautonomic_markers, pt.size = 0, stack = TRUE, flip = TRUE) +
  scale_x_discrete(labels = function(x) sub(" .*", "", x)) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8, hjust = 0.5, angle = 0),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        strip.text.y = element_text(size = 8, face = "bold"))
ggsave(here("plots", "violin_preautonomic_PVN_adult.tiff"), height = 3, width = 5, dpi = 600, compression = "lzw")

# assign PVN/SON NEN identities
adult_PVN_annotations <- tibble(cluster_id = c("C4-293", "C4-295", "C4-297", "C4-298", "C4-312", "C4-313", "C4-314", "C4-315"),
                                NEN_subtype = c("AVP_PVN", "AVP_PVN", "AVP_PVN", "AVP_PVN", "AVP_SON", "OXT", "CRH", "TRH"))

# screen ARC-associated GHRH and dopaminergic populations
adult_ARC_candidates <- adult_NEN_candidates %>%
  filter(features.plot %in% c("GHRH", "TH", "DDC"),
         `Region with highest abdunance score` == "ARC")

# rank neuropeptide expression
adult_ARC_candidates %>%
  select(cluster_id, features.plot, region_threshold_manual,
         `Region with highest abdunance score`, pct.exp, avg.exp) %>%
  arrange(features.plot, desc(avg.exp))

# define ARC neuroendocrine markers
ARC_NEN_markers <- c("GHRH", "TRH", "GAL" , "GSX1", "FOXP2", "DLX1", "DLX2", "TH", "DDC", "SLC6A3", "PRLR", "SLC18A2", "DBH", "NFIB")

# calculate neuroendocrine marker expression across candidate clusters
adult_ARC_expression <- DotPlot(subset(seurat_adult,
                                       cells = colnames(seurat_adult)[seurat_adult$C4_named %in% unique(adult_ARC_candidates$id)]),
                                features = ARC_NEN_markers,
                                group.by = "C4_named",
                                assay = "RNA")$data

# summarise co-marker expression for each cluster
adult_ARC_expression <- adult_ARC_expression %>%
  mutate(cluster_id = sub(" .*", "", as.character(id))) %>%
  select(cluster_id, marker = features.plot, pct.exp, avg.exp) %>%
  pivot_wider(names_from = marker, values_from = c(pct.exp, avg.exp), names_glue = "{marker}_{.value}")

# examine GHRH candidates and co-markers
adult_ARC_expression %>%
  select(cluster_id, starts_with("GHRH_"), starts_with("GSX1_"), starts_with("GHSR_"),
         starts_with("FOXP2_"), starts_with("DLX1_"), starts_with("DLX2_"),
         starts_with("OTP_"), starts_with("SIX6_")) %>%
  arrange(desc(GHRH_avg.exp)) %>%
  print(width = Inf)

# select GHRH candidate populations following neuropeptide, co-marker and anatomical assessment
GHRH_ids <- c("C4-160", "C4-161")

# subset shortlisted populations
GHRH_subset <- subset(seurat_adult, cells = colnames(seurat_adult)[sub(" .*", "", seurat_adult$C4_named) %in% GHRH_ids])
Idents(GHRH_subset) <- "C4_named"

# define GHRH neuroendocrine markers
GHRH_markers <- c("TRH", "GAL", "GSX1", "FOXP2", "DLX1", "DLX2")

# examine neuroendocrine marker expression using stacked violins
VlnPlot(GHRH_subset, features = GHRH_markers, pt.size = 0, stack = TRUE, flip = TRUE)

# assign ARC GHRH NEN identities
adult_GHRH_annotations <- tibble(cluster_id = c("C4-160", "C4-161"),
                                 NEN_subtype = c("GHRH", "GHRH"))

# screen periventricular SST populations
adult_SST_candidates <- adult_NEN_candidates %>%
  filter(features.plot == "SST",
         `Region with highest abdunance score` == "Periventricular")

# rank SST expression
adult_SST_candidates %>%
  select(cluster_id, region_threshold_manual,
         `Region with highest abdunance score`, pct.exp, avg.exp) %>%
  arrange(desc(avg.exp))

# define SST neuroendocrine markers
SST_markers <- c("SST", "PRDM8")

# calculate SST co-marker expression across candidate clusters
adult_SST_expression <- DotPlot(subset(seurat_adult,
                                       cells = colnames(seurat_adult)[seurat_adult$C4_named %in% adult_SST_candidates$id]),
                                features = SST_markers,
                                group.by = "C4_named",
                                assay = "RNA")$data

# summarise co-marker expression for each cluster
adult_SST_expression <- adult_SST_expression %>%
  mutate(cluster_id = sub(" .*", "", as.character(id))) %>%
  select(cluster_id, marker = features.plot, pct.exp, avg.exp) %>%
  pivot_wider(names_from = marker, values_from = c(pct.exp, avg.exp), names_glue = "{marker}_{.value}")

# examine SST candidates and co-markers
adult_SST_expression %>%
  select(cluster_id, starts_with("SST_"), starts_with("PRDM8_")) %>%
  arrange(desc(SST_avg.exp))

# define SST candidates
SST_ids <- c("C4-363")

# subset shortlisted population
SST_subset <- subset(seurat_adult,
                     cells = colnames(seurat_adult)[sub(" .*", "", seurat_adult$C4_named) %in% SST_ids])
Idents(SST_subset) <- "C4_named"

# stacked violin plot
VlnPlot(SST_subset, features = SST_markers, assay = "RNA", pt.size = 0, stack = TRUE, flip = TRUE) 

# screen POA-associated GnRH populations
adult_GnRH_candidates <- adult_NEN_candidates %>%
  filter(features.plot == "GNRH1",
         `Region with highest abdunance score` %in% c("POA", "MPOA", "LPOA"))

# rank GnRH expression
adult_GnRH_candidates %>%
  select(cluster_id, region_threshold_manual,
         `Region with highest abdunance score`, pct.exp, avg.exp) %>%
  arrange(desc(avg.exp))

# define GnRH neuroendocrine markers
GnRH_markers <- c("GNRH1", "KISS1R", "COCH", "CXCL14", "EBF3", "SEMA3C", "SEMA4D", "SOX4", "PDE11A",
                  "FREM1", "MEIS1", "SIX6", "ISL1", "CHAT", "SLC18A3")

# calculate GnRH co-marker expression across candidate clusters
adult_GnRH_expression <- DotPlot(subset(seurat_adult,
                                        cells = colnames(seurat_adult)[seurat_adult$C4_named %in% adult_GnRH_candidates$id]),
                                 features = GnRH_markers,
                                 group.by = "C4_named",
                                 assay = "RNA")$data

# summarise co-marker expression for each cluster
adult_GnRH_expression <- adult_GnRH_expression %>%
  mutate(cluster_id = sub(" .*", "", as.character(id))) %>%
  select(cluster_id, marker = features.plot, pct.exp, avg.exp) %>%
  pivot_wider(names_from = marker, values_from = c(pct.exp, avg.exp), names_glue = "{marker}_{.value}")

# examine GnRH candidates and co-markers
adult_GnRH_expression %>%
  select(cluster_id, starts_with("GNRH1_"), starts_with("KISS1R_"),
         starts_with("COCH_"), starts_with("CXCL14_"), starts_with("EBF3_"),
         starts_with("SEMA3C_"), starts_with("SEMA4D_"), starts_with("SOX4_"),
         starts_with("PDE11A_"), starts_with("FREM1_"), starts_with("MEIS1_"),
         starts_with("SIX6_"), starts_with("ISL1_"),
         starts_with("CHAT_"), starts_with("SLC18A3_")) %>%
  arrange(desc(GNRH1_avg.exp)) %>%
  print(width = Inf)

# combine PVN/SON and ARC NEN annotations
adult_NEN_annotations <- bind_rows(adult_PVN_annotations, adult_GHRH_annotations)
adult_NEN_annotations

# assign NEN subtypes to adult seurat metadata
adult_C4_ids <- sub(" .*", "", seurat_adult$C4_named)
seurat_adult$NEN_subtype <- unname(setNames(adult_NEN_annotations$NEN_subtype, adult_NEN_annotations$cluster_id)[adult_C4_ids])

# classify NENs, non-NEN neurons and non-neuronal cells
seurat_adult$NEN_type <- case_when(!is.na(seurat_adult$NEN_subtype) ~ "NEN", 
                                   seurat_adult$C0_named == "C0-3 Neuron" ~ "non_NEN",
                                   TRUE ~ "non_neuron")

# assign non-NEN and non-neuronal identities
seurat_adult$NEN_subtype[is.na(seurat_adult$NEN_subtype)] <-seurat_adult$NEN_type[is.na(seurat_adult$NEN_subtype)]

# check number of NENs and subtypes
table(seurat_adult$NEN_type)
table(seurat_adult$NEN_subtype)

# verify NEN cluster assignments
table(seurat_adult$C4_named[seurat_adult$NEN_type == "NEN"], seurat_adult$NEN_subtype[seurat_adult$NEN_type == "NEN"])

# define adult NEN subtype order
adult_NEN_order <- c("CRH", "AVP_PVN", "AVP_SON", "OXT", "TRH", "GHRH")

# define colours
adult_NEN_colours <- c("CRH" = "#DD805E", "AVP_PVN" = "#A995C9", "AVP_SON" = "#78A8D1", "OXT" = "#D78FA9", 
                       "TRH" = "#C5A04A", "GHRH" = "#5FA9A2")

# pie chart for proportions of adult NEN subtypes
NEN_subtype_summary <- seurat_adult@meta.data %>%
  filter(NEN_type == "NEN", !is.na(NEN_subtype), !NEN_subtype %in% c("non_NEN", "non_neuron")) %>%
  count(NEN_subtype, name = "cell_count") %>%
  mutate(NEN_subtype = factor(NEN_subtype, levels = adult_NEN_order)) %>%
  arrange(NEN_subtype) %>%
  mutate(percentage = cell_count / sum(cell_count) * 100)

ggplot(NEN_subtype_summary, aes(x = "", y = percentage, fill = NEN_subtype)) +
  geom_col(width = 1) +
  geom_text(aes(x = 1.8, label = paste0(round(percentage, 1), "%")),
            position = position_stack(vjust = 0.5),
            size = 2.5) +
  coord_polar(theta = "y", clip = "off") +
  scale_fill_manual(values = adult_NEN_colours) +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0, "pt"))

ggsave(here("plots", "adult_NEN_pie_proportions.tiff"), height = 2.5, width = 2.5, dpi = 600, compression = "lzw")

# calculate NEN subtype proportions by sex
adult_NEN_sex_summary <- seurat_adult@meta.data %>%
  filter(NEN_type == "NEN") %>%
  count(sex, NEN_subtype, name = "cell_count") %>%
  group_by(sex) %>%
  mutate(percentage = cell_count / sum(cell_count) * 100) %>%
  ungroup()

# set sex order
adult_NEN_sex_summary$sex <- factor(adult_NEN_sex_summary$sex, levels = c("Male", "Female"))

# grouped bar chart of NEN subtype proportions by sex
ggplot(adult_NEN_sex_summary, aes(x = NEN_subtype, y = percentage, fill = sex)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  labs(y = "Proportion of NENs (%)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("Male" = "#78A8D1", "Female" = "#D78FA9")) +
  theme_classic() +
  theme(axis.title = element_text(size = 7),
        axis.text = element_text(size = 7),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        legend.position = "bottom",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.3, "cm"),
        legend.spacing.y = unit(0, "mm"),
        legend.box.spacing = unit(0, "pt"))

ggsave(here("plots", "adult_NEN_sex_grouped.tiff"), height = 2, width = 3.5, dpi = 600, compression = "lzw")

# set adult subtype order in metadata
seurat_adult$NEN_subtype <- factor(seurat_adult$NEN_subtype,
                                   levels = c("CRH", "AVP_PVN", "AVP_SON", "OXT", "TRH", "GHRH", "non_NEN", "non_neuron"))

# subset to NENs only for plotting
NEN_adult <- subset(seurat_adult, subset = NEN_type == "NEN" & NEN_subtype %in% adult_NEN_order)

# dot plot of neuropeptides and co-markers
DotPlot(NEN_adult, assay = "RNA", group.by = "NEN_subtype",
        features = c("DLX1", "FOXP2", "GSX1", "THRB", "NFIB", "FOXP1", "KCNMB4", "TH", "TAC1", "AGTR1", "NR3C1", 
                     "GHRH", "TRH", "OXT", "AVP", "CRH"),
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
        panel.heights = unit(16 * 0.4, "cm"),
        panel.widths = unit(length(adult_NEN_order) * 0.5, "cm"))

ggsave(here("plots", "adult_NEN_dotplot.tiff"), height = 3.1, width = 2.3, dpi = 600, compression = "lzw")

# identify C4 clusters assigned to each anatomical region
adult_region_clusters <- adult_spatial %>%
  filter(grepl("^C4-", `C4 cluster`)) %>%
  dplyr::select(cluster_id = `C4 cluster`, region = `Region with highest abdunance score`)

# define region-specific clusters
POA_ids <- adult_region_clusters %>%
  filter(region %in% c("POA", "MPOA", "LPOA")) %>%
  pull(cluster_id)

ARC_ids <- adult_region_clusters %>%
  filter(region == "ARC") %>%
  pull(cluster_id)

periventricular_ids <- adult_region_clusters %>%
  filter(region == "Periventricular") %>%
  pull(cluster_id)

# POA neurons
POA_subset <- subset(seurat_adult, cells = colnames(seurat_adult)[sub(" .*", "", seurat_adult$C4_named) %in% POA_ids])

# ARC neurons
ARC_subset <- subset(seurat_adult, cells = colnames(seurat_adult)[sub(" .*", "", seurat_adult$C4_named) %in% ARC_ids])

# periventricular neurons
periventricular_subset <- subset(seurat_adult, cells = colnames(seurat_adult)[sub(" .*", "", seurat_adult$C4_named) %in% periventricular_ids])

# POA GnRH expression
VlnPlot(POA_subset, features = "GNRH1", group.by = "C4_named", assay = "RNA", pt.size = 0) +
  scale_x_discrete(labels = function(x) sub(" .*", "", x)) +
  labs(x = NULL, y = "Expression level", title = "POA: GNRH1") +
  theme_classic() +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 7, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        legend.position = "none")
ggsave("violin_POA_GNRH1_adult.tiff", height = 2.5, width = 9.9)

# define GnRH-associated and cholinergic markers
markers_to_plot <- c("GNRH1", "ISL1", "KISS1R", "COCH", "CXCL14", "ACHE", "CHAT", "SLC18A3")

# extract marker expression from C4-194
GnRH_194_markers <- FetchData(GnRH_194, vars = markers_to_plot, layer = "data")

# calculate percentage of C4-194 cells expressing each marker
GnRH_194_barplot <- GnRH_194_markers %>%
  summarise(across(everything(), ~ mean(.x > 0) * 100)) %>%
  pivot_longer(everything(), names_to = "gene", values_to = "percent_expressed") %>%
  mutate(gene = factor(gene, levels = markers_to_plot),
         gene_class = factor(ifelse(gene %in% c("ACHE", "CHAT", "SLC18A3"), "Cholinergic", "GnRH-associated"),
                             levels = c("GnRH-associated", "Cholinergic")))

# plot marker expression across C4-194
ggplot(GnRH_194_barplot, aes(x = gene, y = percent_expressed, fill = gene_class)) +
  geom_col(width = 0.7) +
  facet_grid(~gene_class, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("Cholinergic" = "#78A8D1",
                               "GnRH-associated" = "#A995C9")) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Cells expressing marker (%)") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        legend.position = "none")

ggsave(here("plots", "adult_GnRH_C4-194_markers.tiff"), width = 4.5, height = 2, dpi = 600, compression = "lzw")

# ARC dopaminergic expression
VlnPlot(ARC_subset, features = "TH", group.by = "C4_named", assay = "RNA") +
  scale_x_discrete(labels = function(x) sub(" .*", "", x)) +
  labs(x = NULL, y = "Expression level", title = "ARC: TH") +
  theme_classic() +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 7, hjust = 1, angle = 90),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        legend.position = "none")
ggsave("violin_ARC_TH_adult.tiff", height = 2.5, width = 5.5)

# periventricular somatostatin expression
VlnPlot(periventricular_subset, features = "SST", group.by = "C4_named", assay = "RNA") +
  scale_x_discrete(labels = function(x) sub(" .*", "", x)) +
  labs(x = NULL, y = "Expression level", title = "Periventricular: SST") +
  theme_classic() +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 7, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        legend.position = "none")
ggsave("violin_peri_SST_adult.tiff", height = 2.5, width = 5.5)
