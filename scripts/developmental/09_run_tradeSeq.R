# load packages
library(here)
library(tidyverse)
library(tradeSeq)
library(monocle3)
library(Matrix)
library(BiocParallel)
library(AnnotationDbi)
library(org.Hs.eg.db)

# set working directory
setwd(here())

# get cells from the main connected Monocle trajectory
trade_cells <- colnames(cds_pvn_main)

# extract the raw gene-count matrix for these cells
counts_trade <- counts(cds_pvn_main)[, trade_cells, drop = FALSE]

# number of genes and cells
dim(counts_trade) # 47,429 genes across 19,892 cells

# extract the Monocle pseudotime value for every trajectory cell
pseudotime_vector <- pseudotime(cds_pvn_main)[trade_cells]

# retrieve the nearest graph node for every trajectory cell
trade_cell_nodes <- cell_metadata$closest_node[
  match(trade_cells, cell_metadata$cell)
]

# assign each cell a weight of 1 when its graph node lies on a lineage path
cell_weights <- vapply(endpoint_paths,
                       function(path_nodes) {
                         as.numeric(trade_cell_nodes %in% path_nodes)
                       },
                       numeric(length(trade_cells)))

# label matrix rows and columns
rownames(cell_weights) <- trade_cells
colnames(cell_weights) <- names(endpoint_paths)

# view structure
dim(cell_weights) # 19,892 cells across 3 lineage paths

# names of lineage paths
colnames(cell_weights)

# count how many lineage paths each cell belongs to
table(rowSums(cell_weights)) # 0 - side branches outside selected trajectories, 1 - lineage-specific, 2 - shared by 2 endpoints, 3 - common ancesral trunk shared by all 3

# keep cells belonging to at least one selected lineage
keep_trade_cells <- rowSums(cell_weights) > 0 &
  is.finite(pseudotime_vector)

# retain the same cells, in the same order, across all objects
trade_cells <- trade_cells[keep_trade_cells] # the cell names/barcodes and thier order
counts_trade <- counts_trade[, trade_cells, drop = FALSE] # raw expression counts: genes x cells
cell_weights <- cell_weights[trade_cells, , drop = FALSE] # which lineage path each cell belongs to
pseudotime_vector <- pseudotime_vector[trade_cells] # pseudotime value for each cell

# calculate sequencing-depth offsets using all genes
full_library_sizes <- Matrix::colSums(counts_trade)
trade_offset <- log(pmax(full_library_sizes, 1))

# repeat Monocle pseudotime across the three lineage columns 
trade_pseudotime <- matrix(pseudotime_vector,
                           nrow = length(pseudotime_vector),
                           ncol = ncol(cell_weights),
                           dimnames = list(trade_cells, colnames(cell_weights)))

# check
dim(counts_trade) # 47,429 genes across 13,356 retained cells
dim(trade_pseudotime) 
dim(cell_weights) # 13,3356 cells across 3 lineages

# keep genes detected in at least 20 trajectory cells
keep_genes <- rowSums(counts_trade > 0) >= 20

# see how many genes pass the filter
table(keep_genes)

# retain the expressed genes
counts_trade <- counts_trade[keep_genes, , drop = FALSE]
dim(counts_trade) # 20,718 genes across 13,356 cells

# genes eligible for knot evaluation
pilot_pool <- rownames(counts_trade)

# randomly select the official default of 500 genes
set.seed(123)
pilot_genes_500 <- sample(pilot_pool, size = 500)
pilot_counts_500 <- counts_trade[pilot_genes_500, , drop = FALSE]
dim(pilot_counts_500)

# evaluate 3–7 knots using 500 informative genes
set.seed(123)
k_evaluation_500 <- evaluateK(counts = as.matrix(pilot_counts_500),
                              pseudotime = as.matrix(trade_pseudotime),
                              cellWeights = as.matrix(cell_weights),
                              offset = as.numeric(trade_offset),
                              k = 3:10,
                              nGenes = 500,
                              plot = FALSE,
                              verbose = TRUE,
                              parallel = TRUE,
                              BPPARAM = BiocParallel::MulticoreParam(workers = 4))

# verify successful fits
colSums(is.finite(k_evaluation_500))

# plot the extended comparison
plot_evalutateK_results(k_evaluation_500,k = 3:10) # improvement plateaus after k = 8

# summarise the average AIC improvement at each knot value
mean_aic <- colMeans(k_evaluation_500)
data.frame(nknots = 3:10,
           mean_AIC = mean_aic,
           AIC_improvement = c(NA, -diff(mean_aic))) # small improvement after k = 8

# divide the filtered genes into groups of 500
gene_chunks <- split(rownames(counts_trade), ceiling(seq_along(rownames(counts_trade)) / 500))

# create a folder for the chunk-level results
tradeSeq_output_dir <- "results/tradeSeq_chunks"
dir.create(tradeSeq_output_dir, recursive = TRUE, showWarnings = FALSE)

# fit and test each gene chunk
for (i in seq_along(gene_chunks)) {
  
  output_file <- file.path(tradeSeq_output_dir, sprintf("chunk_%03d.rds", i))
  
  # skip only fully completed and correctly tested chunks
  if (file.exists(output_file)) {
    
    existing_result <- readRDS(output_file)
    
    if (isTRUE(existing_result$complete) && isTRUE(existing_result$association_corrected)) 
    {
      message("Skipping completed chunk ", i, " of ", length(gene_chunks))
      rm(existing_result)
      next
    }
    
    rm(existing_result)
  }
  
  chunk_start_time <- Sys.time()
  message("Starting chunk ", i, " of ", length(gene_chunks), " at ", chunk_start_time)
  
  # extract counts for the genes in this chunk
  chunk_counts_i <- counts_trade[gene_chunks[[i]], , drop = FALSE]
  
  # fit lineage-specific GAMs
  chunk_fit_i <- fitGAM(counts = as.matrix(chunk_counts_i),
                        pseudotime = as.matrix(trade_pseudotime),
                        cellWeights = as.matrix(cell_weights),
                        offset = as.numeric(trade_offset),
                        nknots = 8,
                        verbose = TRUE,
                        parallel = TRUE,
                        BPPARAM = MulticoreParam(workers = 4, progressbar = TRUE))
  
  # test for expression changes along pseudotime
  association_i <- associationTest(chunk_fit_i, l2fc = log2(1.5), inverse = "eigen", contrastType = "start") %>%
    as.data.frame() %>%
    rownames_to_column("gene")
  
  # test whether expression patterns differ between lineages
  pattern_i <- patternTest(chunk_fit_i) %>%
    as.data.frame() %>%
    rownames_to_column("gene")
  
  # test for expression differences between lineage endpoints
  diff_end_i <- diffEndTest(chunk_fit_i) %>%
    as.data.frame() %>%
    rownames_to_column("gene")
  
  # store the results and analysis settings together
  chunk_result <- list(association = association_i, 
                       pattern = pattern_i, 
                       diff_end = diff_end_i,
                       association_corrected = TRUE,
                       complete = TRUE,
                       parameters = list(nknots = 8,
                                         association_l2fc = log2(1.5), 
                                         association_inverse = "eigen",
                                         association_contrastType = "start",
                                         workers = 4))
  
  # save only after every test has completed
  saveRDS(chunk_result, output_file)
  
  chunk_end_time <- Sys.time()
  
  message("Finished and saved chunk ", i, " of ", length(gene_chunks), " | Runtime: ",
          round(as.numeric(difftime(chunk_end_time, chunk_start_time, units = "mins")), 1), " minutes")
  
  # remove large temporary objects
  rm(chunk_counts_i, chunk_fit_i, association_i, pattern_i, diff_end_i, chunk_result)
  gc()
}

# load all 42 standard chunk files
standard_chunk_files <- list.files(tradeSeq_output_dir, pattern = "^chunk_[0-9]{3}\\.rds$", full.names = TRUE)
all_chunk_results <- lapply(standard_chunk_files, readRDS)

# combine pseudotime-association results
association_all <- map_dfr(all_chunk_results, "association") %>%
  transmute(gene, 
            association_waldStat = waldStat, 
            association_df = df,
            association_pvalue = pvalue,
            association_meanLogFC = meanLogFC,
            association_FDR = p.adjust(association_pvalue, method = "BH"))

# combine lineage-pattern results
pattern_all <- map_dfr(all_chunk_results, "pattern") %>%
  transmute(gene, 
            pattern_waldStat = waldStat,
            pattern_df = df,
            pattern_pvalue = pvalue,
            pattern_fcMedian = fcMedian,
            pattern_FDR = p.adjust(pattern_pvalue, method = "BH"))

# combine terminal-endpoint results
diff_end_all <- map_dfr(all_chunk_results, "diff_end") %>%
  transmute(gene, 
            endpoint_waldStat = waldStat, 
            endpoint_df = df,
            endpoint_pvalue = pvalue,
            logFC1_2,
            logFC1_3,
            logFC2_3,
            endpoint_FDR = p.adjust(endpoint_pvalue, method = "BH"),
            max_abs_logFC = pmax(abs(logFC1_2), abs(logFC1_3), abs(logFC2_3), na.rm = TRUE),
            max_abs_logFC = if_else(is.infinite(max_abs_logFC), NA_real_, max_abs_logFC))

# join the three tradeSeq tests
tradeSeq_results <- association_all %>%
  full_join(pattern_all, by = "gene") %>%
  full_join(diff_end_all, by = "gene")

# name the three tradeSeq lineages
lineage_names <- c("OXT/AVP-OXT", "TRH/AVP-CRH", "non-NEN")
colnames(cell_weights) <- lineage_names
colnames(trade_pseudotime) <- lineage_names

# extract metadata in tradeSeq cell order
trade_metadata <- PVN_trajectory_dev@meta.data[rownames(cell_weights), , drop = FALSE] %>%
  as.data.frame()
trade_metadata$cell <- rownames(trade_metadata)

# encode lineage membership
weight_pattern <- apply(1L * (cell_weights > 0), 1, paste0, collapse = "")

# create normalised pseudotime and endpoint matrices
normalised_pseudotime <- matrix(NA_real_, nrow = nrow(trade_pseudotime), ncol = ncol(trade_pseudotime), dimnames = dimnames(trade_pseudotime))
endpoint_membership <- matrix(FALSE, nrow = nrow(trade_pseudotime), ncol = ncol(trade_pseudotime), dimnames = dimnames(trade_pseudotime))

# scale pseudotime and identify the final 10% of each lineage
for (lineage in seq_along(lineage_names)) {
  
  in_lineage <- cell_weights[, lineage] > 0
  lineage_pseudotime <- trade_pseudotime[in_lineage, lineage]
  
  normalised_pseudotime[in_lineage, lineage] <- (lineage_pseudotime - min(lineage_pseudotime)) / 
    (max(lineage_pseudotime) - min(lineage_pseudotime))
  
  endpoint_cutoff <- quantile(lineage_pseudotime, 0.90, na.rm = TRUE)
  
  endpoint_membership[in_lineage, lineage] <- lineage_pseudotime >= endpoint_cutoff
}

# create the final cell-level trajectory metadata
trajectory_state_metadata <- trade_metadata %>%
  mutate(weight_pattern = weight_pattern,
         pseudotime_OXT_AVP_OXT = normalised_pseudotime[, "OXT/AVP-OXT"],
         pseudotime_TRH_AVP_CRH = normalised_pseudotime[, "TRH/AVP-CRH"],
         pseudotime_non_NEN = normalised_pseudotime[, "non-NEN"],
         terminal_OXT_AVP_OXT = endpoint_membership[, "OXT/AVP-OXT"],
         terminal_TRH_AVP_CRH = endpoint_membership[, "TRH/AVP-CRH"],
         terminal_non_NEN = endpoint_membership[, "non-NEN"],
         
         developmental_phase = case_when(trajectory_stage %in% c("Neuroepithelial/RG", "Cycling progenitor") ~ "Early progenitor",
                                         trajectory_stage %in% c("IPC", "Immature neuron") ~ "Intermediate differentiation",
                                         trajectory_stage == "PVN neuron" ~ "PVN neuron", TRUE ~ NA_character_),
         
         trajectory_position = case_when(terminal_OXT_AVP_OXT ~ "Terminal: OXT/AVP-OXT",
                                         terminal_TRH_AVP_CRH ~ "Terminal: TRH/AVP-CRH",
                                         terminal_non_NEN ~ "Terminal: non-NEN",
                                         weight_pattern == "111" ~ "Shared early trunk",
                                         weight_pattern == "100" ~ "OXT/AVP-OXT lineage: non-terminal",
                                         weight_pattern == "011" ~ "Shared TRH/AVP-CRH + non-NEN lineage",
                                         weight_pattern == "010" ~ "TRH/AVP-CRH lineage: non-terminal",
                                         weight_pattern == "001" ~ "non-NEN lineage: non-terminal", 
                                         TRUE ~ NA_character_),
         
         terminal_PVN_group = case_when(trajectory_stage == "PVN neuron" & terminal_OXT_AVP_OXT ~ "OXT/AVP-OXT",
                                        trajectory_stage == "PVN neuron" & terminal_TRH_AVP_CRH ~ "TRH/AVP-CRH",
                                        trajectory_stage == "PVN neuron" & terminal_non_NEN ~ "non-NEN",
                                        TRUE ~ NA_character_))

# check the final assignments
trajectory_state_metadata %>% dplyr::count(developmental_phase, trajectory_position)
table(trajectory_state_metadata$terminal_PVN_group, useNA = "ifany")

# define the terminal PVN groups
terminal_group_order <- c("OXT/AVP-OXT", "TRH/AVP-CRH", "non-NEN")

# identify the cells belonging to each terminal group
terminal_cells <- setNames(lapply(terminal_group_order,
                                  function(group) {
                                    trajectory_state_metadata$cell[which(trajectory_state_metadata$terminal_PVN_group == group)]
                                  }),
                           terminal_group_order)

# calculate the percentage of cells expressing each gene
terminal_detection <- sapply(terminal_cells,
                             function(cells) {
                               rowMeans(counts_trade[, cells, drop = FALSE] > 0) * 100
                             })

# convert the detection matrix into a gene-level table
rownames(terminal_detection) <- rownames(counts_trade)
terminal_detection <- as.data.frame(terminal_detection)
colnames(terminal_detection) <- c("detection_OXT_AVP_OXT", "detection_TRH_AVP_CRH", "detection_non_NEN")
terminal_detection <- terminal_detection %>%
  rownames_to_column("gene") %>%
  mutate(max_terminal_detection = pmax(detection_OXT_AVP_OXT, detection_TRH_AVP_CRH, detection_non_NEN),
         min_terminal_detection = pmin(detection_OXT_AVP_OXT, detection_TRH_AVP_CRH, detection_non_NEN))

# add terminal detection to the master tradeSeq table
tradeSeq_results <- tradeSeq_results %>%
  left_join(terminal_detection, by = "gene")

# retrieve current NCBI gene types
ncbi_gene_types <- mapIds(org.Hs.eg.db,
                          keys = unique(tradeSeq_results$gene),
                          column = "GENETYPE",
                          keytype = "SYMBOL",
                          multiVals = "first")

# download the curated human transcription factor list
human_tfs <- readLines("https://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt") %>%
  trimws() %>%
  unique()

# filter for protein-coding TFs and non-TFs only
tradeSeq_results <- tradeSeq_results %>%
  mutate(ncbi_gene_type = unname(ncbi_gene_types[gene]),
         gene_class = case_when(ncbi_gene_type == "protein-coding" & gene %in% human_tfs ~ "TF",
                                ncbi_gene_type == "protein-coding" & !gene %in% human_tfs ~ "non-TF",
                                TRUE ~ "Excluded"))

# filter out ribosomal, mitochondiral and haemoglobin genes for plotting
tradeSeq_results <- tradeSeq_results %>%
  mutate(technical_gene = case_when(
    
    # cytosolic ribosomal genes
    str_detect(gene, "^(RPL([0-9]+[A-Z]{0,2}|P[0-9]+)|RPS([0-9]+[A-Z]{0,2}|A))$") ~ TRUE,
    
    # mitochondrially encoded genes
    str_detect(gene, "^MT-") ~ TRUE,
    
    # haemoglobin genes
    gene %in% c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ") ~ TRUE,
    
    TRUE ~ FALSE))

# classify terminal expression programmes
tradeSeq_results <- tradeSeq_results %>%
  mutate(terminal_program = case_when(
    
    # branch-specific programmes
    endpoint_FDR < 0.001 &
      max_terminal_detection >= 5 &
      logFC1_2 >= 2 &
      logFC1_3 >= 2 ~
      "Strong OXT/AVP-OXT-enriched",
    
    endpoint_FDR < 0.001 &
      max_terminal_detection >= 5 &
      logFC1_2 <= -2 &
      logFC2_3 >= 2 ~
      "Strong TRH/AVP-CRH-enriched",
    
    endpoint_FDR < 0.001 &
      max_terminal_detection >= 5 &
      logFC1_3 <= -2 &
      logFC2_3 <= -2 ~
      "Strong non-NEN-enriched",
    
    # programmes shared by the two NEN branches
    endpoint_FDR < 0.001 &
      max_terminal_detection >= 5 &
      abs(logFC1_2) < 1 &
      logFC1_3 >= 1 &
      logFC2_3 >= 1 ~
      "Shared NEN-enriched",
    
    # programmes enriched in non-NEN cells
    endpoint_FDR < 0.001 &
      max_terminal_detection >= 5 &
      abs(logFC1_2) < 1 &
      logFC1_3 <= -1 &
      logFC2_3 <= -1 ~
      "non-NEN-enriched",
    
    # genes expressed similarly across all three endpoints
    min_terminal_detection >= 5 &
      max_abs_logFC < 1 ~
      "Shared across all terminals",
    
    # significant divergence without one clean pattern
    endpoint_FDR < 0.001 &
      max_terminal_detection >= 5 &
      max_abs_logFC >= 2 ~
      "Other terminal divergence",
    
    TRUE ~
      "Not terminally prioritised"))

# protein-coding transcription factors
tf_results <- tradeSeq_results %>%
  filter(gene_class == "TF", !technical_gene, !str_detect(gene, fixed("///")))

# protein-coding non-transcription factors
non_tf_results <- tradeSeq_results %>%
  filter(gene_class == "non-TF", !technical_gene, !str_detect(gene, fixed("///")))

# check the classifications
bind_rows(tf_results %>%
            dplyr::count(terminal_program) %>%
            mutate(gene_class = "TF"),
          
          non_tf_results %>%
            dplyr::count(terminal_program) %>%
            mutate(gene_class = "non-TF")) %>%
  arrange(gene_class, desc(n))

# retain terminal PVN cells
terminal_metadata <- trajectory_state_metadata %>%
  filter(!is.na(terminal_PVN_group), cell %in% colnames(counts_trade))

# extract terminal cells as columns
terminal_counts <- counts_trade[, terminal_metadata$cell, drop = FALSE]

# log-normalise counts
terminal_logcounts <- terminal_counts %*%
  Diagonal(x = 10000 / colSums(terminal_counts))
terminal_logcounts@x <- log1p(terminal_logcounts@x)

# restore gene and cell names after matrix normalisation
dimnames(terminal_logcounts) <- dimnames(terminal_counts)

# calculate mean expression in each endpoint
terminal_mean_expression <- sapply(terminal_groups,
                                   function(group) {
                                     group_cells <- terminal_metadata$cell[terminal_metadata$terminal_PVN_group == group]
                                     rowMeans(terminal_logcounts[, group_cells, drop = FALSE])
                                   })
colnames(terminal_mean_expression) <-
  terminal_groups

# terminal programmes displayed in the TF heatmap
tf_program_levels <- c("OXT/AVP-OXT-enriched", "TRH/AVP-CRH-enriched", "Shared NEN-enriched", "non-NEN-enriched", "Shared across all terminals")

# retain all TFs assigned to a clear terminal programme
tf_divergent_genes <- tf_results %>%
  filter(terminal_program %in% c("Strong OXT/AVP-OXT-enriched", "Strong TRH/AVP-CRH-enriched", "Strong non-NEN-enriched",
                                 "Shared NEN-enriched", "non-NEN-enriched"))

# select broadly detected TFs with minimal terminal differences
tf_shared_genes <- tf_results %>%
  filter(terminal_program == "Shared across all terminals",
         max_abs_logFC < 0.5, min_terminal_detection >= 20) %>%
  arrange(max_abs_logFC,  desc(min_terminal_detection)) %>%
  slice_head(n = 10)