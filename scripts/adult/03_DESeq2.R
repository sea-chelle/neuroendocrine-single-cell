library(Seurat)
library(DESeq2)
library(tidyverse)

# create donor-population metadata
pb_meta <- neuron_adult@meta.data %>%
  dplyr::count(Donor_ID, sex, NEN_subtype, name = "n_cells") %>%
  mutate(pb_id = paste0("PB", row_number()))

# assign pseudobulk IDs to nuclei
neuron_adult$pb_id <- pb_meta$pb_id[match(paste(neuron_adult$Donor_ID, neuron_adult$NEN_subtype),
                                          paste(pb_meta$Donor_ID, pb_meta$NEN_subtype))]

# aggregate raw counts by donor and population
pb_counts <- AggregateExpression(neuron_adult,
                                 assays = "RNA",
                                 group.by = "pb_id",
                                 return.seurat = FALSE)$RNA

# group pseudobulk samples into NEN and non_NEN
pb_meta <- pb_meta %>%
  mutate(NEN_type = ifelse(NEN_subtype == "non_NEN", "non_NEN", "NEN"))

# combine NEN subtype counts within each donor
NEN_pb_counts <- sapply(split(pb_meta$pb_id, paste(pb_meta$Donor_ID, pb_meta$NEN_type, sep = ".")),
                        function(ids) rowSums(pb_counts[, ids, drop = FALSE]))

# create combined pseudobulk metadata
NEN_pb_meta <- pb_meta %>%
  group_by(Donor_ID, sex, NEN_type) %>%
  summarise(n_cells = sum(n_cells), .groups = "drop") %>%
  mutate(pb_id = paste(Donor_ID, NEN_type, sep = ".")) %>%
  as.data.frame()

# align metadata with count matrix
rownames(NEN_pb_meta) <- NEN_pb_meta$pb_id
NEN_pb_meta <- NEN_pb_meta[colnames(NEN_pb_counts), , drop = FALSE]

# verify alignment
stopifnot(identical(rownames(NEN_pb_meta), colnames(NEN_pb_counts)))

# retain donors with at least 20 nuclei in both populations
eligible_donors <- NEN_pb_meta %>%
  filter(n_cells >= 20) %>%
  dplyr::count(Donor_ID) %>%
  filter(n == 2) %>%
  pull(Donor_ID)

# filter and prepare pseudobulk metadata
NEN_pb_meta_filtered <- NEN_pb_meta %>%
  filter(Donor_ID %in% eligible_donors) %>%
  mutate(Donor_ID = factor(Donor_ID),
         sex = factor(sex, levels = c("Female", "Male")),
         NEN_type = factor(NEN_type, levels = c("non_NEN", "NEN"))) %>%
  as.data.frame()

rownames(NEN_pb_meta_filtered) <- NEN_pb_meta_filtered$pb_id

# select corresponding pseudobulk counts
NEN_pb_counts_filtered <- NEN_pb_counts[, rownames(NEN_pb_meta_filtered), drop = FALSE]

# verify alignment
stopifnot(identical(colnames(NEN_pb_counts_filtered), rownames(NEN_pb_meta_filtered)))

# define function for DESeq2 analysis
run_pseudobulk_DESeq2 <- function(counts, meta, design, contrast, min_samples = 3) {
  
  # align metadata and counts
  meta <- droplevels(meta[colnames(counts), , drop = FALSE])
  stopifnot(identical(rownames(meta), colnames(counts)))
  
  # create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = round(counts), colData = meta, design = design)
  
  # filter lowly expressed genes
  keep <- rowSums(DESeq2::counts(dds) >= 10) >= min_samples
  dds <- dds[keep, ]
  
  # run differential expression analysis
  dds <- DESeq(dds)
  
  # extract results
  res <- results(dds, contrast = contrast) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    arrange(padj) %>%
    mutate(diffexpressed = case_when(padj < 0.05 & log2FoldChange > 0.5 ~ "UP", 
                                     padj < 0.05 & log2FoldChange < -0.5 ~ "DOWN",
                                     TRUE ~ "NO"))
  
# define NEN subtypes
NEN_subtypes <- setdiff(unique(pb_meta$NEN_subtype), "non_NEN")

# function to compare each NEN subtype against remaining NENs
run_NEN_subtype_DESeq2 <- function(subtype) {

  # select NEN pseudobulk samples
  meta <- pb_meta %>%
    filter(NEN_subtype != "non_NEN") %>%
    mutate(Group = ifelse(NEN_subtype == subtype, subtype, "Other_NEN"))

  # aggregate counts by donor and group
  group_ids <- split(meta$pb_id, paste(meta$Donor_ID, meta$Group, sep = "."))

  counts <- sapply(group_ids, function(ids) {
    rowSums(pb_counts[, ids, drop = FALSE])
  })

  # create pseudobulk metadata
  group_meta <- meta %>%
    group_by(Donor_ID, sex, Group) %>%
    summarise(n_cells = sum(n_cells), .groups = "drop") %>%
    mutate(pb_id = paste(Donor_ID, Group, sep = ".")) %>%
    as.data.frame()

  rownames(group_meta) <- group_meta$pb_id

  # retain donors with at least 20 nuclei in both groups
  eligible_donors <- group_meta %>%
    filter(n_cells >= 20) %>%
    dplyr::count(Donor_ID) %>%
    filter(n == 2) %>%
    pull(Donor_ID)

  group_meta <- group_meta %>%
    filter(Donor_ID %in% eligible_donors) %>%
    mutate(Donor_ID = factor(Donor_ID),
           Group = factor(Group, levels = c("Other_NEN", subtype))) %>%
    as.data.frame()

  rownames(group_meta) <- group_meta$pb_id

  # select corresponding counts
  counts <- counts[, rownames(group_meta), drop = FALSE]

  # run DESeq2
  run_pseudobulk_DESeq2(counts = counts,
                        meta = group_meta,
                        design = ~ Donor_ID + Group,
                        contrast = c("Group", subtype, "Other_NEN"),
                        min_samples = 3)
}
  
  return(list(dds = dds, results = res))
}

# define DESeq2 comparisons
NEN_comparisons <- list(Overall = list(subset = rep(TRUE, nrow(NEN_pb_meta_filtered)),
                                       design = ~ Donor_ID + NEN_type,
                                       contrast = c("NEN_type", "NEN", "non_NEN"),
                                       min_samples = 4),
                        
                        Male_vs_Female = list(subset = NEN_pb_meta_filtered$NEN_type == "NEN",
                                              design = ~ sex,
                                              contrast = c("sex", "Male", "Female"),
                                              min_samples = 2),
                        
                        Male_NEN = list(subset = NEN_pb_meta_filtered$sex == "Male",
                                        design = ~ Donor_ID + NEN_type,
                                        contrast = c("NEN_type", "NEN", "non_NEN"),
                                        min_samples = 3),
                        
                        Female_NEN = list(subset = NEN_pb_meta_filtered$sex == "Female",
                                          design = ~ Donor_ID + NEN_type,
                                          contrast = c("NEN_type", "NEN", "non_NEN"),
                                          min_samples = 2))

# run all DESeq2 comparisons
NEN_analyses <- lapply(names(NEN_comparisons), function(name) {
  
  config <- NEN_comparisons[[name]]
  
  # subset metadata and counts
  meta <- NEN_pb_meta_filtered[config$subset, , drop = FALSE]
  counts <- NEN_pb_counts_filtered[, meta$pb_id, drop = FALSE]
  
  # run DESeq2
  run_pseudobulk_DESeq2(counts = counts,
                        meta = meta,
                        design = config$design,
                        contrast = config$contrast,
                        min_samples = config$min_samples)
})

names(NEN_analyses) <- names(NEN_comparisons)

# extract DEG results
NEN_DE_results <- lapply(NEN_analyses, `[[`, "results")

# summarise DEGs across comparisons
NEN_DE_summary <- imap_dfr(NEN_DE_results, function(res, name) {
  data.frame(comparison = name,
             UP = sum(res$diffexpressed == "UP"),
             DOWN = sum(res$diffexpressed == "DOWN"),
             NO = sum(res$diffexpressed == "NO"),
             total_DEGs = sum(res$diffexpressed != "NO"))
})

NEN_DE_summary

# summarise eligible paired donors for each NEN subtype sex-DEG comparison
map_dfr(NEN_subtypes, function(subtype) {
  
  pb_meta %>%
    filter(NEN_subtype != "non_NEN") %>%
    mutate(Group = ifelse(NEN_subtype == subtype, subtype, "Other_NEN")) %>%
    group_by(Donor_ID, sex, Group) %>%
    summarise(n_cells = sum(n_cells), .groups = "drop") %>%
    filter(n_cells >= 20) %>%
    dplyr::count(Donor_ID, sex) %>%
    filter(n == 2) %>%
    summarise(subtype = subtype,
              n_donors = n(),
              n_male = sum(sex == "Male"),
              n_female = sum(sex == "Female"))
})

# define NEN subtypes
NEN_subtypes <- setdiff(unique(pb_meta$NEN_subtype), "non_NEN")

# function to compare each NEN subtype against remaining NENs
run_NEN_subtype_DESeq2 <- function(subtype) {
  
  # select NEN pseudobulk samples
  meta <- pb_meta %>%
    filter(NEN_subtype != "non_NEN") %>%
    mutate(Group = ifelse(NEN_subtype == subtype, subtype, "Other_NEN"))
  
  # aggregate counts by donor and group
  group_ids <- split(meta$pb_id, paste(meta$Donor_ID, meta$Group, sep = "."))
  
  counts <- sapply(group_ids, function(ids) {
    rowSums(pb_counts[, ids, drop = FALSE])
  })
  
  # create pseudobulk metadata
  group_meta <- meta %>%
    group_by(Donor_ID, sex, Group) %>%
    summarise(n_cells = sum(n_cells), .groups = "drop") %>%
    mutate(pb_id = paste(Donor_ID, Group, sep = ".")) %>%
    as.data.frame()
  
  rownames(group_meta) <- group_meta$pb_id
  
  # retain donors with at least 20 nuclei in both groups
  eligible_donors <- group_meta %>%
    filter(n_cells >= 20) %>%
    dplyr::count(Donor_ID) %>%
    filter(n == 2) %>%
    pull(Donor_ID)
  
  group_meta <- group_meta %>%
    filter(Donor_ID %in% eligible_donors) %>%
    mutate(Donor_ID = factor(Donor_ID),
           Group = factor(Group, levels = c("Other_NEN", subtype))) %>%
    as.data.frame()
  
  rownames(group_meta) <- group_meta$pb_id
  
  # select corresponding counts
  counts <- counts[, rownames(group_meta), drop = FALSE]
  
  # run DESeq2
  run_pseudobulk_DESeq2(counts = counts,
                        meta = group_meta,
                        design = ~ Donor_ID + Group,
                        contrast = c("Group", subtype, "Other_NEN"),
                        min_samples = 3)
}

# run all NEN subtype comparisons
NEN_subtype_analyses <- lapply(NEN_subtypes, run_NEN_subtype_DESeq2)
names(NEN_subtype_analyses) <- NEN_subtypes

# extract DEG results
NEN_subtype_DE_results <- lapply(NEN_subtype_analyses, `[[`, "results")

# summarise DEGs across subtype comparisons
NEN_subtype_DE_summary <- imap_dfr(NEN_subtype_DE_results, function(res, name) {
  data.frame(subtype = name,
             UP = sum(res$diffexpressed == "UP"),
             DOWN = sum(res$diffexpressed == "DOWN"),
             NO = sum(res$diffexpressed == "NO"),
             total_DEGs = sum(res$diffexpressed != "NO"))
})

NEN_subtype_DE_summary

# extract top 20 upregulated genes for each NEN subtype
NEN_subtype_top_genes <- purrr::imap_dfr(NEN_subtype_DE_results, function(res, subtype) {
  res %>%
    filter(diffexpressed == "UP") %>%
    arrange(desc(log2FoldChange)) %>%
    slice_head(n = 20) %>%
    mutate(subtype = subtype)
})

# display top genes by subtype
NEN_subtype_top_genes %>%
  dplyr::select(subtype, gene, log2FoldChange, padj)

# summarise developmental NEN subtype representation
nuclei_dev@meta.data %>%
  filter(NEN_subtype != "non_NEN") %>%
  dplyr::count(NEN_subtype, donor, name = "n_cells") %>%
  group_by(NEN_subtype) %>%
  summarise(n_donors = n_distinct(donor),
            n_donors_20 = sum(n_cells >= 20),
            total_cells = sum(n_cells),
            min_cells = min(n_cells),
            max_cells = max(n_cells),
            .groups = "drop") %>%
  arrange(desc(total_cells))

# create developmental donor-population metadata
dev_pb_meta <- nuclei_dev@meta.data %>%
  dplyr::count(donor, NEN_subtype, name = "n_cells") %>%
  mutate(pb_id = paste0("Dev_PB", row_number()))

# assign pseudobulk IDs to nuclei
nuclei_dev$pb_id <- dev_pb_meta$pb_id[match(paste(nuclei_dev$donor, nuclei_dev$NEN_subtype),
                                            paste(dev_pb_meta$donor, dev_pb_meta$NEN_subtype))]

# aggregate raw counts by donor and population
dev_pb_counts <- AggregateExpression(nuclei_dev,
                                     assays = "RNA",
                                     group.by = "pb_id",
                                     return.seurat = FALSE)$RNA

# restore pseudobulk IDs
colnames(dev_pb_counts) <- gsub("-", "_", colnames(dev_pb_counts))

# verify pseudobulk IDs
stopifnot(all(dev_pb_meta$pb_id %in% colnames(dev_pb_counts)))

# define developmental vers us adult NEN comparisons
dev_adult_comparisons <- list(Overall = list(dev = c("AVP-CRH", "AVP-OXT", "OXT", "TRH", "GHRH"),
                                             adult = c("AVP_PVN", "AVP_SON", "CRH", "OXT", "TRH", "GHRH")),
                              GHRH = list(dev = "GHRH", adult = "GHRH"),
                              TRH = list(dev = "TRH", adult = "TRH"),
                              AVP_CRH_vs_CRH = list(dev = "AVP-CRH", adult = "CRH"),
                              AVP_CRH_vs_AVP_PVN = list(dev = "AVP-CRH", adult = "AVP_PVN"),
                              AVP_OXT_vs_OXT = list(dev = "AVP-OXT", adult = "OXT"),
                              AVP_OXT_vs_AVP_PVN = list(dev = "AVP-OXT", adult = "AVP_PVN"))

# function to compare developmental and adult NEN populations
run_dev_adult_DESeq2 <- function(dev_subtypes, adult_subtypes) {
  
  # select developmental pseudobulk samples
  dev_meta <- dev_pb_meta %>%
    filter(as.character(NEN_subtype) %in% dev_subtypes) %>%
    group_by(donor) %>%
    summarise(n_cells = sum(n_cells),
              ids = list(pb_id),
              .groups = "drop") %>%
    filter(n_cells >= 20)
  
  # select adult pseudobulk samples
  adult_meta <- pb_meta %>%
    filter(NEN_subtype %in% adult_subtypes) %>%
    group_by(Donor_ID) %>%
    summarise(n_cells = sum(n_cells),
              ids = list(pb_id),
              .groups = "drop") %>%
    filter(n_cells >= 20)
  
  # check biological replication
  if (nrow(dev_meta) < 2 || nrow(adult_meta) < 2) {
    stop("Insufficient donors for comparison")
  }
  
  # aggregate developmental counts by donor
  dev_counts <- sapply(dev_meta$ids, function(ids) {
    rowSums(dev_pb_counts[, ids, drop = FALSE])
  })
  
  colnames(dev_counts) <- paste0("Dev_", dev_meta$donor)
  
  # aggregate adult counts by donor
  adult_counts <- sapply(adult_meta$ids, function(ids) {
    rowSums(pb_counts[, ids, drop = FALSE])
  })
  
  colnames(adult_counts) <- paste0("Adult_", adult_meta$Donor_ID)
  
  # retain genes present in both datasets
  common_genes <- intersect(rownames(dev_counts), rownames(adult_counts))
  
  counts <- cbind(dev_counts[common_genes, , drop = FALSE],
                  adult_counts[common_genes, , drop = FALSE])
  
  # create combined pseudobulk metadata
  meta <- data.frame(pb_id = colnames(counts),
                     Donor_ID = c(as.character(dev_meta$donor), as.character(adult_meta$Donor_ID)),
                     Stage = factor(c(rep("Developmental", nrow(dev_meta)), rep("Adult", nrow(adult_meta))), 
                                    levels = c("Developmental", "Adult")),
                     n_cells = c(dev_meta$n_cells, adult_meta$n_cells))
  
  rownames(meta) <- meta$pb_id
  
  # verify alignment
  stopifnot(identical(colnames(counts), rownames(meta)))
  
  # run DESeq2
  analysis <- run_pseudobulk_DESeq2(counts = counts,
                                    meta = meta,
                                    design = ~ Stage,
                                    contrast = c("Stage", "Adult", "Developmental"),
                                    min_samples = 2)
  
  # retain sample metadata
  analysis$meta <- meta
  
  return(analysis)
}

# run developmental versus adult DESeq2 comparisons
dev_adult_analyses <- lapply(names(dev_adult_comparisons), function(name) {
  
  config <- dev_adult_comparisons[[name]]
  
  run_dev_adult_DESeq2(dev_subtypes = config$dev,
                       adult_subtypes = config$adult)
})

names(dev_adult_analyses) <- names(dev_adult_comparisons)

# extract DEG results
dev_adult_DE_results <- lapply(dev_adult_analyses, `[[`, "results")

# summarise DEGs across comparisons
dev_adult_DE_summary <- purrr::imap_dfr(dev_adult_analyses, function(analysis, name) {
  
  res <- analysis$results
  
  data.frame(comparison = name,
             n_dev = sum(analysis$meta$Stage == "Developmental"),
             n_adult = sum(analysis$meta$Stage == "Adult"),
             UP = sum(res$diffexpressed == "UP"),
             DOWN = sum(res$diffexpressed == "DOWN"),
             NO = sum(res$diffexpressed == "NO"),
             total_DEGs = sum(res$diffexpressed != "NO"))
})

dev_adult_DE_summary

# extract top 20 DEGs for each developmental-adult comparison
dev_adult_top_genes <- purrr::imap_dfr(dev_adult_DE_results, function(res, name) {
  bind_rows(res %>%
      filter(diffexpressed == "UP") %>%
      arrange(desc(log2FoldChange)) %>%
      slice_head(n = 20) %>%
      mutate(direction = "Adult"),
    res %>%
      filter(diffexpressed == "DOWN") %>%
      arrange(log2FoldChange) %>%
      slice_head(n = 20) %>%
      mutate(direction = "Developmental")) %>%
    mutate(comparison = name)
})

# examine global NEN DEGs
dev_adult_top_genes %>%
  filter(comparison == "Overall") %>%
  dplyr::select(comparison, direction, gene, log2FoldChange, padj)