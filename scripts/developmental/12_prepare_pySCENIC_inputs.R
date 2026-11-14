# load required packages
library(here)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(Matrix)

# vreate a dedicated output directory for SCENIC
scenic_output_dir <- file.path(here("results", "pySCENIC", "dev"))
dir.create(scenic_output_dir, recursive = TRUE, showWarnings = FALSE)

# define all immature neurons
immature_cells <- colnames(neuron_dev)[as.character(neuron_dev$cell_type) %in% "Immature"]

# define mature PVN and ARC neurons
pvn_cells <- colnames(nuclei_dev)[nuclei_dev$NEN_region %in% "PVN"]
arc_cells <- colnames(nuclei_dev)[nuclei_dev$NEN_region %in% "ARC"]

# all mature neurons of interest
mature_cells <- unique(c(pvn_cells, arc_cells))

# complete population in which regulon activity will be scored
scenic_scoring_cells <- unique(c(progenitor_cells, immature_cells, mature_cells))

# extract mature neuronal annotations
mature_annotations <- nuclei_dev@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(cell, 
            mature_region = as.character(NEN_region),
            NEN_status_original = as.character(NEN_type),
            NEN_subtype = as.character(NEN_subtype))

# prepare metadata in scoring matrix order
scenic_metadata <- seurat_dev@meta.data[scenic_scoring_cells, , drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  left_join(mature_annotations, by = "cell") %>%
  mutate(developmental_stage = case_when(cell %in% mature_cells ~ "Mature neuron",
                                         cell %in% immature_cells ~ "Immature neuron",
                                         as.character(seurat_clusters) == "22" ~ "Neuroepithelial/RG",
                                         as.character(seurat_clusters) %in% c("25", "32", "39", "40") ~ "Cycling progenitor",
                                         as.character(seurat_clusters) == "13" ~ "IPC",
                                          TRUE ~ NA_character_),
         maturity = case_when(developmental_stage == "Mature neuron" ~ "Mature",
                              developmental_stage %in% c("Neuroepithelial/RG", "Cycling progenitor", "IPC", "Immature neuron") ~ "Developing",
                              TRUE ~ NA_character_),
         NEN_status = case_when(maturity == "Mature" & NEN_status_original == "NEN" ~ "NEN",
                                maturity == "Mature" & NEN_status_original == "non_NEN" ~ "non-NEN",
                                TRUE ~ NA_character_),
         scenic_group = case_when(maturity == "Developing" ~ developmental_stage,
                                  mature_region == "PVN" & NEN_subtype == "OXT" ~ "PVN: OXT",
                                  mature_region == "PVN" & NEN_subtype == "AVP-OXT" ~ "PVN: AVP-OXT",
                                  mature_region == "PVN" & NEN_subtype == "TRH" ~ "PVN: TRH",
                                  mature_region == "PVN" & NEN_subtype == "AVP-CRH" ~ "PVN: AVP-CRH",
                                  mature_region == "PVN" & NEN_status == "non-NEN" ~ "PVN: non-NEN",
                                  mature_region == "ARC" & NEN_subtype == "GHRH" ~ "ARC: GHRH",
                                  mature_region == "ARC" & NEN_status == "non-NEN" ~ "ARC: non-NEN",
                                  TRUE ~ NA_character_))

# inspect group composition
table(scenic_metadata$scenic_group, useNA = "ifany")

# sample up to 2,000 cells from each biological group
set.seed(123)
scenic_inference_metadata <- scenic_metadata %>% group_by(scenic_group) %>%
  group_modify(~ slice_sample(.x, n = min(2000, nrow(.x)))) %>%
  ungroup()

scenic_inference_cells <- scenic_inference_metadata$cell

# record inference membership in the complete metadata
scenic_metadata <- scenic_metadata %>%
  mutate(used_for_inference = cell %in% scenic_inference_cells)

# check number of cells in each group used for inference
scenic_inference_metadata %>%
  dplyr::count(scenic_group, name = "inference_cells") %>%
  arrange(desc(inference_cells))

# extract raw RNA counts
scenic_raw_counts <- LayerData(seurat_dev, assay = "RNA", layer = "counts")

# inference matrix in exactly the saved cell order
scenic_inference_counts <- scenic_raw_counts[, scenic_inference_cells, drop = FALSE]

# retain genes detected in at least 0.5% of inference cells
minimum_detection_fraction <- 0.005
minimum_cells <- ceiling(ncol(scenic_inference_counts) * minimum_detection_fraction)

# require the equivalent of at least three counts per minimum cell number
minimum_total_counts <- 3 * minimum_cells
c(minimum_cells = minimum_cells, minimum_total_counts = minimum_total_counts)

# calculate gene-level expression summaries
gene_filter_table <- tibble(gene = rownames(scenic_inference_counts),
                            detected_cells = rowSums(scenic_inference_counts > 0),
                            total_counts = rowSums(scenic_inference_counts)) %>%
  mutate(technical_gene = str_detect(gene, "^(MT-|RPL[0-9A-Z]+$|RPS[0-9A-Z]+$)") | 
           gene %in% c("MALAT1", "XIST", "HBA1","HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1",  "HBZ") |
           str_detect(gene, fixed("///")),
         duplicated_gene = duplicated(gene),
         retained = detected_cells >= minimum_cells & total_counts >= minimum_total_counts & !technical_gene & !duplicated_gene)

# inspect filtering outcome
gene_filter_table %>%
  dplyr::count(retained)

# retain genes in their original matrix order
scenic_genes <- gene_filter_table %>%
  filter(retained) %>%
  pull(gene)

length(scenic_genes)

# balanced matrix used for GRNBoost2 and motif pruning
scenic_inference_counts <- scenic_inference_counts[scenic_genes, , drop = FALSE]

# complete matrix used for AUCell scoring
scenic_scoring_counts <- scenic_raw_counts[, scenic_scoring_cells, drop = FALSE]

# check final dimensions
c(inference_genes = nrow(scenic_inference_counts), 
  inference_cells = ncol(scenic_inference_counts),
  scoring_genes = nrow(scenic_scoring_counts),
  scoring_cells = ncol(scenic_scoring_counts))

# restore exact metadata order 
scenic_metadata <- scenic_metadata %>%
  slice(match(scenic_scoring_cells, cell))

scenic_inference_metadata <- scenic_metadata %>%
  slice(match(scenic_inference_cells, cell))

# save metadata and filtering records
write_csv(scenic_metadata, file.path(scenic_output_dir, "pyscenic_scoring_metadata.csv"))
write_csv(scenic_inference_metadata, file.path(scenic_output_dir, "pyscenic_inference_metadata.csv" ))
write_csv(gene_filter_table, file.path(scenic_output_dir, "pyscenic_gene_filter.csv"))

# write and save filter paramters
filter_parameters <- tibble(parameter = c("random_seed", "maximum_cells_per_group", "minimum_detection_fraction", 
                                          "minimum_detected_cells", "minimum_total_counts", "retained_genes", "inference_cells", "scoring_cells"),
                            value = c(123, 2000, minimum_detection_fraction, minimum_cells, minimum_total_counts,
                                      length(scenic_genes), length(scenic_inference_cells), length(scenic_scoring_cells)))


write_csv(filter_parameters, file.path(scenic_output_dir, "pyscenic_filter_parameters.csv"))

# save ordered gene and cell identifiers
writeLines(scenic_genes, file.path(scenic_output_dir, "pyscenic_inference_genes.txt"))
writeLines(rownames(scenic_scoring_counts), file.path(scenic_output_dir, "pyscenic_scoring_genes.txt"))
writeLines(scenic_inference_cells, file.path(scenic_output_dir, "pyscenic_inference_cells.txt"))
writeLines(scenic_scoring_cells, file.path(scenic_output_dir, "pyscenic_scoring_cells.txt"))

# export sparse matrices for Python 
writeMM(scenic_inference_counts, file.path(scenic_output_dir, "pyscenic_inference_counts.mtx"))
writeMM(scenic_scoring_counts, file.path(scenic_output_dir, "pyscenic_scoring_counts.mtx"))

# check
stopifnot(identical(colnames(scenic_inference_counts), scenic_inference_cells),
          identical(colnames(scenic_scoring_counts), scenic_scoring_cells),
          identical(scenic_metadata$cell, scenic_scoring_cells),
          identical(scenic_inference_metadata$cell, scenic_inference_cells),
          !anyNA(scenic_metadata$scenic_group),
          !anyDuplicated(rownames(scenic_scoring_counts)))