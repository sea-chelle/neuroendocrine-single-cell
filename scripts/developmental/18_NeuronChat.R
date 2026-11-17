# load requied packages
library(here)
library(Seurat)
library(NeuronChat)
library(tidyverse)
library(purrr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# set working directory 
setwd(here("plots"))

# use all annotated hypothalamic neurons
neuronchat_adult <- subset(neuron_adult, subset = !spatial_region %in% c("Thalamus", "Thalamaus", "Ventricular", "NA"))

# retain neuroendocrine subtype annotations and use nuclei annotations for other neurons
neuronchat_adult$neuronchat_group <- ifelse(neuronchat_adult$NEN_type == "NEN",
                                            as.character(neuronchat_adult$NEN_subtype),
                                            as.character(neuronchat_adult$spatial_region))

# distinguish neuroendocrine subtypes by their hypothalamic nuclei
neuronchat_adult$neuronchat_group <- recode(neuronchat_adult$neuronchat_group, 
                                            "CRH" = "PVN: CRH",
                                            "AVP_PVN" = "PVN: AVP",
                                            "AVP_SON" = "SON: AVP",
                                            "OXT" = "PVN: OXT",
                                            "TRH" = "PVN: TRH",
                                            "GHRH" = "ARC: GHRH",
                                            "PVN" = "PVN: non-NEN",
                                            "ARC" = "ARC: non-NEN")

# inspect population sizes
table(neuronchat_adult$neuronchat_group, useNA = "ifany")

# extract normalised counts
neuronchat_data <- GetAssayData(neuronchat_adult, assay = "RNA", layer = "data")

# cell group labels for NeuronChat
neuronchat_labels <- neuronchat_adult$neuronchat_group

# make sure labels are named by cell barcode
names(neuronchat_labels) <- colnames(neuronchat_adult)

# check cell type labels are in the same order as neuronchat data barcodes
identical(colnames(neuronchat_data), names(neuronchat_labels))

# create neuronchat object
neuronchat <- createNeuronChat(object = neuronchat_data, group.by = neuronchat_labels, DB = "human")

# run neuronchat
neuronchat <- run_NeuronChat(neuronchat)

# extract all inferred neuronal communication interactions
neuronchat_links <- imap_dfr(neuronchat@net, function(mat, interaction) {
  as.data.frame(as.table(mat)) %>%
    setNames(c("sender", "receiver", "strength")) %>%
    mutate(interaction = interaction)
}) %>%
  mutate(sender = as.character(sender),
         receiver = as.character(receiver)) %>%
  filter(strength > 0)

# count nonzero communication links for each interaction
neuronchat_network_summary <- tibble(interaction = names(neuronchat@net),
                                     n_links = sapply(neuronchat@net, function(x) sum(x > 0, na.rm = TRUE))) %>%
  arrange(desc(n_links))

# check number of interactions with inferred communication
neuronchat_network_summary %>%
  summarise(total_interactions = n(),
            active_interactions = sum(n_links > 0),
            total_links = sum(n_links))

# summarise inferred communication between every pair of neuronal populations
neuronchat_pair_summary <- neuronchat_links %>%
  group_by(sender, receiver) %>%
  summarise(n_interactions = n_distinct(interaction),
            total_strength = sum(strength),
            mean_strength = mean(strength),
            .groups = "drop") %>%
  arrange(desc(total_strength))

# inspect strongest predicted sender-receiver pairs
head(neuronchat_pair_summary, 30) %>% 
  print(n = Inf)

# summarise total incoming communication for each neuronal population
neuronchat_incoming <- neuronchat_links %>%
  group_by(receiver) %>%
  summarise(incoming_links = n(),
            incoming_interactions = n_distinct(interaction),
            incoming_strength = sum(strength),
            incoming_partners = n_distinct(sender),
            .groups = "drop") %>%
  rename(population = receiver)

# summarise total outgoing communication for each neuronal population
neuronchat_outgoing <- neuronchat_links %>%
  group_by(sender) %>%
  summarise(outgoing_links = n(),
            outgoing_interactions = n_distinct(interaction),
          outgoing_strength = sum(strength),
          outgoing_partners = n_distinct(receiver),
          .groups = "drop") %>%
  rename(population = sender)

# combine incoming and outgoing communication summaries
neuronchat_communication_summary <- neuronchat_incoming %>%
  full_join(neuronchat_outgoing, by = "population")

# inspect communication profiles across all populations
print(neuronchat_communication_summary, n = Inf, width = Inf)

# define neuroendocrine populations of interest
neuronchat_NEN_groups <- c("PVN: AVP-CRH", "PVN: AVP-OXT", "PVN: OXT", "PVN: TRH", "ARC: GHRH")

# compare communication profiles across the five NEN subtypes
neuronchat_communication_summary %>%
  filter(population %in% neuronchat_NEN_groups) %>%
  print(width = Inf)

# define PVN neuroendocrine subtypes
neuronchat_PVN_NEN <- c("PVN: AVP-CRH", "PVN: AVP-OXT", "PVN: OXT", "PVN: TRH")

# combine PVN NEN subtypes for the global network visualisation
neuronchat_combined_links <- neuronchat_links %>%
  mutate(sender = if_else(sender %in% neuronchat_PVN_NEN, "PVN NEN", sender),
    receiver = if_else(receiver %in% neuronchat_PVN_NEN, "PVN NEN", receiver)) %>%
  group_by(sender, receiver) %>%
  summarise(total_strength = sum(strength), n_links = n(), .groups = "drop")

# define the neuroendocrine populations of interest
neuronchat_combined_NEN <- c("PVN NEN", "ARC: GHRH")

# rank other neuronal populations by their communication with NENs
neuronchat_top_partners <- neuronchat_combined_links %>%
  filter(sender %in% neuronchat_combined_NEN |
           receiver %in% neuronchat_combined_NEN) %>%
  filter(!(sender %in% neuronchat_combined_NEN &
             receiver %in% neuronchat_combined_NEN)) %>%
  mutate(partner = if_else(sender %in% neuronchat_combined_NEN,
                           receiver, sender)) %>%
  group_by(partner) %>%
  summarise(total_strength = sum(total_strength), .groups = "drop") %>%
  arrange(desc(total_strength)) %>%
  slice_head(n = 6) %>%
  pull(partner)

# inspect selected populations
neuronchat_top_partners

# define plotting populations
neuronchat_focus_groups <- c(neuronchat_combined_NEN, neuronchat_top_partners)

# retain communication involving NENs and selected partners
neuronchat_focus_links <- neuronchat_combined_links %>%
  filter(sender %in% neuronchat_focus_groups,
         receiver %in% neuronchat_focus_groups,
         sender %in% neuronchat_combined_NEN | receiver %in% neuronchat_combined_NEN)

# create complete sender-receiver matrix
neuronchat_focus_matrix <- matrix(0, 
                                  nrow = length(neuronchat_focus_groups),
                                  ncol = length(neuronchat_focus_groups),
                                  dimnames = list(neuronchat_focus_groups, neuronchat_focus_groups))

# add communication strengths
neuronchat_focus_matrix[cbind(match(neuronchat_focus_links$sender, neuronchat_focus_groups), 
                              match(neuronchat_focus_links$receiver, neuronchat_focus_groups))] <- neuronchat_focus_links$total_strength

# define colours
neuronchat_combined_colours <- c("PVN NEN" = "#A995C9",
                                 "ARC: GHRH" = "#5FA9A2",
                                 "ARC: non-NEN" = "#78A8D1",
                                 "Tuberal" = "#D78FA9",
                                 "TMN" = "#C5A04A",
                                 "LH" = "#DD805E",
                                 "VMH" = "#B85F5A",
                                 "SCN" = "#829B73")

# plot NEN-associated neuronal communication network
pdf("NeuronChat_NEN_circle.pdf", width = 2.5, height = 3, useDingbats = FALSE)

# reset circular plotting parameters
circos.clear()
circos.par(gap.degree = 5, canvas.xlim = c(-1.1, 1.1), canvas.ylim = c(-1.1, 1.1))

# plot communication network
chordDiagram(x = neuronchat_focus_matrix,
             order = neuronchat_focus_groups,
             grid.col = neuronchat_combined_colours,
             transparency = 0.7,
             directional = 1,
             direction.type = "arrows",
             link.arr.type = "big.arrow",
             link.border = NA,
             annotationTrack = "grid",
             annotationTrackHeight = mm_h(1),
             preAllocateTracks = list(track.height = mm_h(3)))

# add neuronal population labels outside the circle
circos.trackPlotRegion(track.index = 1,  bg.border = NA,
                       panel.fun = function(x, y) {
                         sector_name <- get.cell.meta.data("sector.index")
                         
                         # shorten / wrap labels
                         sector_name <- sub("^ARC: GHRH$", "GHRH", sector_name)
                         sector_name <- sub("^ARC: non-NEN$", "ARC\nnon-\nNEN", sector_name)
                         sector_name <- sub("^PVN NEN$", "PVN\nNEN", sector_name) 
                         xlim <- get.cell.meta.data("xlim")
                         ylim <- get.cell.meta.data("ylim")
                         
                         circos.text(x = mean(xlim),
                                     y = ylim[1] + 1.2,
                                     labels = sector_name,
                                     facing = "clockwise",
                                     niceFacing = TRUE,
                                     adj = c(0.5, 0.5),
                                     cex = 0.6)
                         })

circos.clear()
dev.off()

# summarise total communication strength for every sender-receiver pair
neuronchat_pair_strength <- neuronchat_links %>%
  group_by(sender, receiver) %>%
  summarise(total_strength = sum(strength), .groups = "drop")

# define all neuronal populations in the original network
neuronchat_all_groups <- rownames(neuronchat@net[[1]])

# construct the complete communication matrix
neuronchat_strength_matrix <- matrix(0,
                                     nrow = length(neuronchat_all_groups),
                                     ncol = length(neuronchat_all_groups),
                                     dimnames = list(neuronchat_all_groups, neuronchat_all_groups))

# fill the matrix with summed communication strengths
neuronchat_strength_matrix[cbind(match(neuronchat_pair_strength$sender, neuronchat_all_groups),
                                 match(neuronchat_pair_strength$receiver, neuronchat_all_groups))] <- neuronchat_pair_strength$total_strength

# extract incoming communication: all neurons → NENs
neuronchat_incoming_matrix <- neuronchat_strength_matrix[neuronchat_all_groups, neuronchat_NEN_groups, drop = FALSE]

# extract outgoing communication: NENs → all neurons
neuronchat_outgoing_matrix <- t(neuronchat_strength_matrix[neuronchat_NEN_groups, neuronchat_all_groups, drop = FALSE])

# define a shared colour scale across incoming and outgoing communication
neuronchat_max_strength <- max(neuronchat_incoming_matrix, neuronchat_outgoing_matrix)
neuronchat_heatmap_colours <- colorRamp2(c(0, neuronchat_max_strength / 2, neuronchat_max_strength), 
                                         c("white", "#A995C9", "#7656AB"))

# create shorter labels for plotting
short_group_labels <- c("PVN: AVP-CRH" = "\nAVP-\nCRH", "PVN: AVP-OXT" = "\nAVP-\nOXT", "PVN: OXT" = "OXT",
                        "PVN: TRH" = "TRH", "PVN: non-NEN" = "PVN\nnon-NEN", "ARC: GHRH" = "GHRH", "ARC: non-NEN" = "ARC\nnon-NEN")

# shorten population labels in both heatmap matrices
rownames(neuronchat_incoming_matrix) <- recode(rownames(neuronchat_incoming_matrix), !!!short_group_labels)
colnames(neuronchat_incoming_matrix) <- recode(colnames(neuronchat_incoming_matrix), !!!short_group_labels)
rownames(neuronchat_outgoing_matrix) <- recode(rownames(neuronchat_outgoing_matrix), !!!short_group_labels)
colnames(neuronchat_outgoing_matrix) <- recode(colnames(neuronchat_outgoing_matrix), !!!short_group_labels)

# incoming communication heatmap
neuronchat_incoming_heatmap <- Heatmap(neuronchat_incoming_matrix,
                                       name = "Communication\nstrength",
                                       col = neuronchat_heatmap_colours,
                                       cluster_rows = FALSE,
                                       cluster_columns = FALSE,
                                       row_names_gp = gpar(fontsize = 7),
                                       row_names_centered = TRUE,
                                       column_names_gp = gpar(fontsize = 7),
                                       column_names_rot = 0,
                                       column_title = "Incoming:\nneurons → NENs",
                                       column_title_gp = gpar(fontsize = 8),
                                       column_names_centered = TRUE,
                                       heatmap_legend_param = list(title_gp = gpar(fontsize = 7),
                                                                   labels_gp = gpar(fontsize = 7),
                                                                   grid_width = unit(3, "mm"),
                                                                   grid_height = unit(5, "mm"),
                                                                   at = c(0, 0.2, 0.4)))

# outgoing communication heatmap
neuronchat_outgoing_heatmap <- Heatmap(neuronchat_outgoing_matrix,
                                       name = "Communication\nstrength",
                                       col = neuronchat_heatmap_colours,
                                       cluster_rows = FALSE,
                                       cluster_columns = FALSE,
                                       row_names_gp = gpar(fontsize = 7),
                                       row_names_centered = TRUE,
                                       column_names_gp = gpar(fontsize = 7),
                                       column_names_rot = 0,
                                       column_title = "Outgoing:\nNENs → neurons",
                                       column_title_gp = gpar(fontsize = 8),
                                       column_names_centered = TRUE,
                                       heatmap_legend_param = list(title_gp = gpar(fontsize = 7),
                                                                   labels_gp = gpar(fontsize = 7),
                                                                   grid_width = unit(3, "mm"),
                                                                   grid_height = unit(5, "mm"),
                                                                   at = c(0, 0.2, 0.4)))

# display the heatmaps together
pdf("neuronchat_heatmap.pdf", height = 4, width = 4.27, useDingbats = FALSE)
draw(neuronchat_incoming_heatmap + neuronchat_outgoing_heatmap, merge_legend = TRUE)
dev.off()

# summarise incoming communication by NEN subtype and molecular interaction
neuronchat_incoming_interactions <- neuronchat_links %>%
  filter(receiver %in% neuronchat_NEN_groups) %>%
  group_by(NEN = receiver, interaction) %>%
  summarise(total_strength = sum(strength), n_partners = n_distinct(sender), .groups = "drop") %>%
  mutate(direction = "Incoming")

# summarise outgoing communication by NEN subtype and molecular interaction
neuronchat_outgoing_interactions <- neuronchat_links %>%
  filter(sender %in% neuronchat_NEN_groups) %>%
  group_by(NEN = sender, interaction) %>%
  summarise(total_strength = sum(strength), n_partners = n_distinct(receiver), .groups = "drop") %>%
  mutate(direction = "Outgoing")

# combine incoming and outgoing interactions
neuronchat_interaction_summary <- bind_rows( neuronchat_incoming_interactions, neuronchat_outgoing_interactions)

# inspect the results
head(neuronchat_interaction_summary, 20)

# order interactions according to their total inferred communication strength
neuronchat_interaction_order <- neuronchat_interaction_summary %>%
  group_by(interaction) %>%
  summarise(total_strength = sum(total_strength), .groups = "drop") %>%
  arrange(desc(total_strength)) %>%
  pull(interaction)

# define NEN subtype order
neuronchat_NEN_order <- c("PVN: AVP-CRH", "PVN: AVP-OXT", "PVN: OXT",  "PVN: TRH", "ARC: GHRH")

# prepare factors for plotting
neuronchat_interaction_plot <- neuronchat_interaction_summary %>%
  mutate(NEN = factor(NEN, levels = neuronchat_NEN_order),
         interaction = factor(interaction, levels = rev(neuronchat_interaction_order)),
         direction = factor(direction, levels = c("Incoming", "Outgoing")))

# plot molecular interactions associated with NEN communication
ggplot(neuronchat_interaction_plot, aes(x = NEN, y = interaction)) +
  geom_point(aes(size = n_partners, colour = total_strength)) +
  facet_wrap(~direction, ncol = 2) +
  scale_x_discrete(labels = function(x) {
    x <- sub("^(PVN|ARC): ", "", x)
    sub("^(AVP)-(CRH|OXT)$", "\\1\n\\2", x)
  }) +
  scale_size_continuous(name = "Communication\npartners", range = c(1, 4) ) +
  scale_colour_gradient(low = "#E3DAF0", high = "#7656AB", name = "Total communication\nstrength") +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.position = "right",
        legend.key.size = unit(3, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.box.spacing = unit(0, "pt"),
        panel.spacing = unit(0.5, "lines"),
        panel.heights = unit(22 * 0.4, "cm"))
ggsave("neurochat_dotplot.tiff", width = 6.27, height = 4, dpi = 600, compression = "lzw")

# identify neuronal populations with predicted NRXN1-NLGN1 signalling to AVP-CRH
neuronchat_links %>%
  filter(interaction == "NRXN1_NLGN1", receiver == "PVN: AVP-CRH") %>%
  arrange(desc(strength)) %>%
  select(sender, receiver, strength)

neuronchat_links %>%
  filter(interaction == "NRXN1_NLGN1", sender == "PVN: AVP-CRH") %>%
  arrange(desc(strength)) %>%
  select(sender, receiver, strength)
