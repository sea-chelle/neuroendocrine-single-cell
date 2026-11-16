# load required packages
library(here)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(Matrix)

# create a dedicated output directory for SCENIC
scenic_output_dir <- file.path(here("results", "pySCENIC", "adult"))
dir.create(scenic_output_dir, recursive = TRUE, showWarnings = FALSE)

# get C4 clusters annotations from the PVN and ARC
adult_region_clusters <- adult_spatial %>%
  filter(grepl("^C4-", `C4 cluster`)) %>%
  mutate(cluster_id = `C4 cluster`, region = region_threshold_manual)

# map spatial region assignments onto adult neurons
adult_cluster_ids <- sub(" .*", "", as.character(neuron_adult$C4_named))
neuron_adult$spatial_region <- adult_region_clusters$region[match(adult_cluster_ids, adult_region_clusters$cluster_id)]

# define PVN and ARC non-NEN populations
pvn_nonNEN_cells <- colnames(neuron_adult)[neuron_adult$spatial_region %in% "PVN" & neuron_adult$NEN_type %in% "non_NEN"]
arc_nonNEN_cells <- colnames(neuron_adult)[neuron_adult$spatial_region %in% "ARC" & neuron_adult$NEN_type %in% "non_NEN"]
son_nonNEN_cells <- colnames(neuron_adult)[neuron_adult$spatial_region %in% "SON" & neuron_adult$NEN_type %in% "non_NEN"]

# define all adult NEN populations
nen_cells <- colnames(neuron_adult)[neuron_adult$NEN_type %in% "NEN"]

# combine cells for regulon scoring
scenic_scoring_cells <- unique(c(pvn_nonNEN_cells, arc_nonNEN_cells, son_nonNEN_cells, nen_cells))

# extract metadata in exact scoring matrix order
scenic_metadata <- neuron_adult@meta.data[scenic_scoring_cells, , drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  transmute(cell, 
            sex = as.character(sex),
            donor = as.character(Donor_ID),
            spatial_region = as.character(spatial_region),
            NEN_status = as.character(NEN_type),
            NEN_subtype = as.character(NEN_subtype),
            scenic_group = case_when(NEN_type == "NEN" ~ as.character(NEN_subtype),
                                     spatial_region == "PVN" & NEN_type == "non_NEN" ~ "PVN: non-NEN",
                                     spatial_region == "ARC" & NEN_type == "non_NEN" ~ "ARC: non-NEN",
                                     spatial_region == "SON" & NEN_type == "non_NEN" ~ "SON: non-NEN",
                                     TRUE ~ NA_character_),
            sex_NEN_group = paste(sex, NEN_status, sep = "_"))

# check final groups
table(scenic_metadata$scenic_group, useNA = "ifany")

# check sex distribution across NEN vs non_NEN
table(scenic_metadata$NEN_status, scenic_metadata$sex)

# sample up to 2,000 cells per biological group, retaining all female cells
set.seed(123)

scenic_inference_metadata <- scenic_metadata %>%
  group_by(scenic_group) %>%
  group_modify(~ {
    female <- .x %>% filter(sex == "Female")
    male <- .x %>% filter(sex == "Male")
    remaining <- max(0, 2000 - nrow(female))
    bind_rows(female, slice_sample(male, n = min(remaining, nrow(male))))
  }) %>%
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
scenic_raw_counts <- LayerData(neuron_adult, assay = "RNA", layer = "counts")

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
           gene %in% c("MALAT1", "HBA1","HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1",  "HBZ") |
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

# export sparse matrices
writeMM(scenic_inference_counts, file.path(scenic_output_dir, "pyscenic_inference_counts.mtx"))
writeMM(scenic_scoring_counts, file.path(scenic_output_dir, "pyscenic_scoring_counts.mtx"))

# check
stopifnot(identical(rownames(scenic_inference_counts), scenic_genes),
          identical(colnames(scenic_inference_counts), scenic_inference_cells),
          identical(colnames(scenic_scoring_counts), scenic_scoring_cells),
          identical(scenic_metadata$cell, scenic_scoring_cells),
          identical(scenic_inference_metadata$cell, scenic_inference_cells),
          !anyNA(scenic_metadata$scenic_group),
          !anyDuplicated(rownames(scenic_scoring_counts)))