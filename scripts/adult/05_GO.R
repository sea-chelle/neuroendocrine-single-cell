# load required packages
library(here)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# combine all DESeq2 comparisons with unique names
GO_DE_results <- c(setNames(NEN_DE_results, paste0("Adult_", names(NEN_DE_results))),
                   setNames(NEN_subtype_DE_results, paste0("Subtype_", names(NEN_subtype_DE_results))),
                   setNames(dev_adult_DE_results, paste0("Development_", names(dev_adult_DE_results))))

names(GO_DE_results)

# define genes to exclude from GO enrichment
exclude_genes <- function(genes) {
  
  # mitochondrial-encoded genes
  mitochondrial <- grepl("^MT-", genes)
  
  # haemoglobin genes
  haemoglobin <- genes %in% c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBZ")
  
  # retain remaining genes
  genes[!mitochondrial & !haemoglobin]
}

# function to perform GO Biological Process enrichment
run_GO_enrichment <- function(res, direction) {
  
  # exclude genes with unavailable adjusted p-values
  res <- res %>%
    filter(!is.na(padj), !is.na(gene)) %>%
    filter(gene %in% exclude_genes(gene))
  
  # define background and significant genes
  background <- unique(res$gene)
  
  genes <- res %>%
    filter(diffexpressed == direction) %>%
    pull(gene) %>%
    unique()
  
  # skip comparisons with insufficient DEGs
  if (length(genes) < 5) return(NULL)
  
  # perform GO enrichment
  GO <- enrichGO(gene = genes, 
                 OrgDb = org.Hs.eg.db,
                 keyType = "SYMBOL",
                 ont = "BP",
                 universe = background,
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2,
                 minGSSize = 10,
                 maxGSSize = 500,
                 readable = FALSE)
  
  return(GO)
}

# perform GO enrichment for UP and DOWN genes
GO_results <- lapply(GO_DE_results, function(res) {
  
  list(UP = run_GO_enrichment(res, "UP"), DOWN = run_GO_enrichment(res, "DOWN"))
  
})

# check comparisons
names(GO_results)

# summarise number of significant GO terms across all comparisons
imap_dfr(GO_results, function(results, comparison) {
  
  tibble(comparison = comparison,
         n_UP = nrow(as.data.frame(results$UP)),
         n_DOWN = nrow(as.data.frame(results$DOWN)),
         n_total = n_UP + n_DOWN)
})