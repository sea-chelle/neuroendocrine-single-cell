# subset object for NENs only
NEN_dev <- subset(nuclei_dev, subset = NEN_type == "NEN")

# set idents to NEN subtype
Idents(NEN_dev) <- "NEN_subtype"

# on vs all subtype differentially expressed markers
NEN_markers <- FindAllMarkers(NEN_dev, assay = "RNA", test.use = "MAST",
                              only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25,
                              latent.vars = "nFeature_RNA")

# get top 20 most differentially expressed genes for each one
NEN_top_markers <- NEN_markers %>%
  
  # retain reasonably robust subtype-enriched markers
  filter( p_val_adj < 0.05, avg_log2FC > 0.5, pct.1 >= 0.20, pct.1 - pct.2 >= 0.10) %>%
  
  # select top 20 markers per NEN subtype
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(cluster, desc(avg_log2FC))

# view
NEN_top_markers %>%
  select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
  print(n = Inf)

celltype_week_prop <- NEN_dev@meta.data %>%
  count(gestational.week, NEN_subtype, name = "cell_count") %>%
  group_by(gestational.week) %>%
  mutate(
    total_cells = sum(cell_count),
    proportion = cell_count / total_cells,
    percentage = proportion * 100
  ) %>%
  ungroup() %>%
  mutate( gestational.week = forcats::fct_reorder(gestational.week, readr::parse_number(gestational.week)))

ggplot(
  celltype_week_prop,
  aes(x = gestational.week, y = proportion, fill = NEN_subtype)
) +
  geom_col(width = 0.8) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Gestational week",
    y = "Proportion of cells",
    fill = "Cell type"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8)
  )
