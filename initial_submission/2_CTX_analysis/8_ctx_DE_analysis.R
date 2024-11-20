library(Seurat)
library(DESeq2)
library(tidyverse)
library(dplyr)
library(magrittr)

xenium_ctx <- readRDS("/xenium_Line61_final/ctx/xenium_CTX_metadata_pSyn.rds")

################################################################################
################################################################################
################################################################################

# Filtering main objects to only have cells for specific comparisons of interest

Idents(xenium_ctx) <- "has_pSyn"
neg_cells_ctx <- WhichCells(xenium_ctx, idents = "FALSE")
pos_cells_ctx <- WhichCells(xenium_ctx, idents = "TRUE")


# retaining nTg and pSyn+ Tg cells

xenium_ctx_pSyn_pos <- xenium_ctx[, !colnames(xenium_ctx) %in% neg_cells_ctx]
saveRDS(xenium_ctx_pSyn_pos, file = "/xenium_Line61_final/ctx/xenium_CTX_pSyn_pos_for_DESeq.rds")


# retaining nTg and pSyn- Tg cells

xenium_ctx_pSyn_neg <- xenium_ctx[, !colnames(xenium_ctx) %in% pos_cells_ctx]
saveRDS(xenium_ctx_pSyn_neg, file = "/xenium_Line61_final/ctx/xenium_CTX_pSyn_neg_for_DESeq.rds")


# retaining Tg cells (for pSyn+ vs. pSyn-)

Idents(xenium_ctx) <- "genotype"
tg_cells_ctx <- WhichCells(xenium_ctx, idents = "Tg")
# using subset_opt function found here: https://github.com/alikhuseynov/add-on_R/blob/develop/R/subset_obj_seurat_v2.R 
# because Seurat is bugged and won't let you subset if there aren't some cells from every fov in the subset object
xenium_ctx_tg <- subset_opt(xenium_ctx, cells = tg_cells_ctx) 
saveRDS(xenium_ctx_tg, file = "/xenium_Line61_final/ctx/xenium_CTX_tg_only_for_DESeq.rds")


################################################################################
################################################################################
################################################################################

# General stuff--reading in lists of genes, functions, etc. 

# genes to keep (when filtering DE results for visualization)

genes_to_keep <- read.csv2("/xenium_Line61_final/xenium_upload1/general/genes_to_keep_DE.csv", header = F)
genes_to_keep <- genes_to_keep$V1

# defining selected_celltypes to use later (see fxn and later uses)

ctx_ExNs <- c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", "L5 ET ExN", "L5 NP ExN", "L6 IT ExN", "L6 CT ExN", "L6b ExN")

# Function to perform DESeq2 analysis for comparisons b/w nTg and Tg cells on every celltype defined in selected_celltypes --
# returns all genes, regardless of pval

run_DESeq_nTg_Tg <- function(matrix_list, selected_celltypes) {
  results_list <- list()
  
  for (matrix_name in selected_celltypes) {
    if (matrix_name %in% names(matrix_list)) {
      matrix_data <- matrix_list[[matrix_name]]
      
      # Set coldata
      colData <- data.frame(samples = colnames(matrix_data))
      colData <- colData %>% 
        mutate(genotype = ifelse(grepl("nTg", samples), "nTg", "Tg")) %>%
        column_to_rownames("samples")
      
      # Make DESeq dataset and run DESeq
      dds <- DESeqDataSetFromMatrix(countData = matrix_data,
                                    colData = colData,
                                    design = ~ genotype)
      dds <- DESeq(dds)
      
      # Extract DE results
      result_name <- paste0("result_", matrix_name)
      result <- as.data.frame(results(dds, name = "genotype_Tg_vs_nTg"))
      
      # Store results in the results_list
      results_list[[result_name]] <- result
    } else {
      warning(paste("Matrix", matrix_name, "not found in the provided list. Skipping."))
    }
  }
  
  # Combine all results into a single dataframe
  combined_results <- do.call(rbind, results_list)
  return(combined_results)
}


# Function to perform DESeq2 analysis for comparisons w/in only Tg cells (comparing pSyn+ to pSyn-) on every celltype present --
# returns all genes, regardless of pval

run_DESeq_Tg_only <- function(matrix_list, selected_names) {
  results_list <- list()
  
  for (matrix_name in selected_names) {
    if (matrix_name %in% names(matrix_list)) {
      matrix_data <- matrix_list[[matrix_name]]
      
      # Set coldata (w/ pSyn status)
      colData <- data.frame(samples = colnames(matrix_data))
      colData <- colData %>% 
        mutate(pSyn_status = ifelse(grepl("TRUE", samples), "TRUE", "FALSE")) %>%
        column_to_rownames("samples")
      
      # Make DESeq dataset and run DESeq (using pSyn status as design factor)
      dds <- DESeqDataSetFromMatrix(countData = matrix_data,
                                    colData = colData,
                                    design = ~ pSyn_status)
      dds <- DESeq(dds)
      
      result_name <- paste0("result_", matrix_name)
      result <- as.data.frame(results(dds, name = "pSyn_status_TRUE_vs_FALSE"))
      
      results_list[[result_name]] <- result
    } else {
      warning(paste("Matrix", matrix_name, "not found in the provided list. Skipping."))
    }
  }
  
  # Combine all results into a single dataframe
  combined_results <- do.call(rbind, results_list)
  return(combined_results)
}


################################################################################
################################################################################
################################################################################

# running diffexp and saving results

################################################################################

# 1. nTg vs. pSyn+ Tg

xenium_ctx_pSyn_pos <- readRDS("/xenium_Line61_final/ctx/xenium_CTX_pSyn_pos_for_DESeq.rds")

cts1 <- AggregateExpression(xenium_ctx_pSyn_pos,
                            group.by = c("celltype", "orig.ident"),
                            assays = "Xenium",
                            slot = "counts",
                            return.seurat = F)
# extract matrix
cts1 <- as.data.frame(cts1[["Xenium"]])

# transpose
cts1.t <- as.data.frame(t(cts1))

# get values where to split
splitrows1 <- gsub('_.*', '', rownames(cts1.t))

# split dataframe
cts1_split <- split.data.frame(cts1.t, 
                               f = factor(splitrows1))

# modify list by changing rownames of elements and transposing
cts1_split_modified <- lapply(cts1_split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
})


res_ctx_nTg_vs_pSyn_pos <- run_DESeq_nTg_Tg(cts1_split_modified, ctx_ExNs)
res_ctx_nTg_vs_pSyn_pos$genes <- rownames(res_ctx_nTg_vs_pSyn_pos)
res_ctx_nTg_vs_pSyn_pos$genes <- gsub("^result_", "", res_ctx_nTg_vs_pSyn_pos$genes)
res_ctx_nTg_vs_pSyn_pos <- separate(res_ctx_nTg_vs_pSyn_pos, genes, into = c("celltype", "gene"), sep = "\\.")

res_ctx_nTg_vs_pSyn_pos.df <- as.data.frame(res_ctx_nTg_vs_pSyn_pos)
saveRDS(res_ctx_nTg_vs_pSyn_pos.df, file = "/xenium_Line61_final/ctx/outs/DE_ctx_nTg_vs_pSyn_pos_ALL.rds")


################################################################################

# 2. nTg vs. pSyn- Tg

xenium_ctx_pSyn_neg <- readRDS("/xenium_Line61_final/ctx/xenium_CTX_pSyn_neg_for_DESeq.rds")

cts2 <- AggregateExpression(xenium_ctx_pSyn_neg,
                            group.by = c("celltype", "orig.ident"),
                            assays = "Xenium",
                            slot = "counts",
                            return.seurat = F)
cts2 <- as.data.frame(cts2[["Xenium"]])
cts2.t <- as.data.frame(t(cts2))
splitrows2 <- gsub('_.*', '', rownames(cts2.t))
cts2_split <- split.data.frame(cts2.t, 
                               f = factor(splitrows2))
cts2_split_modified <- lapply(cts2_split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
})


res_ctx_nTg_vs_pSyn_neg <- run_DESeq_nTg_Tg(cts2_split_modified, ctx_ExNs)
res_ctx_nTg_vs_pSyn_neg$genes <- rownames(res_ctx_nTg_vs_pSyn_neg)
res_ctx_nTg_vs_pSyn_neg$genes <- gsub("^result_", "", res_ctx_nTg_vs_pSyn_neg$genes)
res_ctx_nTg_vs_pSyn_neg <- separate(res_ctx_nTg_vs_pSyn_neg, genes, into = c("celltype", "gene"), sep = "\\.")

res_ctx_nTg_vs_pSyn_neg.df <- as.data.frame(res_ctx_nTg_vs_pSyn_neg)
saveRDS(res_ctx_nTg_vs_pSyn_neg.df, file = "/xenium_Line61_final/ctx/outs/DE_ctx_nTg_vs_pSyn_neg_ALL.rds")

################################################################################

# 3. Tg pSyn+ vs. pSyn-

xenium_tg_only <- readRDS("/xenium_Line61_final/ctx/xenium_CTX_tg_only_for_DESeq.rds")

cts3 <- AggregateExpression(xenium_tg_only,
                            group.by = c("celltype", "orig.ident", "has_pSyn"),
                            assays = "Xenium",
                            slot = "counts",
                            return.seurat = F)
cts3 <- as.data.frame(cts3[["Xenium"]])
cts3.t <- as.data.frame(t(cts3))
splitrows3 <- gsub('_.*', '', rownames(cts3.t))
cts3_split <- split.data.frame(cts3.t, 
                               f = factor(splitrows3))
cts3_split_modified <- lapply(cts3_split, function(x){
  rownames(x) <- gsub('.*_(.*_.*)', '\\1', rownames(x))
  t(x)
})


res_Tg_pos_vs_neg <- run_DESeq_Tg_only(cts3_split_modified, ctx_ExNs)
res_Tg_pos_vs_neg$genes <- rownames(res_Tg_pos_vs_neg)
res_Tg_pos_vs_neg$genes <- gsub("^result_", "", res_Tg_pos_vs_neg$genes)
res_Tg_pos_vs_neg <- separate(res_Tg_pos_vs_neg, genes, into = c("celltype", "gene"), sep = "\\.")

res_Tg_pos_vs_neg.df <- as.data.frame(res_Tg_pos_vs_neg)
saveRDS(res_Tg_pos_vs_neg.df, file = "/xenium_Line61_final/ctx/outs/DE_ctx_Tg_only_ALL.rds")


################################################################################
################################################################################
################################################################################

# combining / manipulating dataframes for deciding which genes to use for visualization

res_ctx_nTg_vs_pSyn_pos.df <- readRDS("/xenium_Line61_final/ctx/outs/DE_ctx_nTg_vs_pSyn_pos_ALL.rds")
res_ctx_nTg_vs_pSyn_neg.df <- readRDS("/xenium_Line61_final/ctx/outs/DE_ctx_nTg_vs_pSyn_neg_ALL.rds")


nTg_vs_pSyn_pos_SIGS.df <- res_ctx_nTg_vs_pSyn_pos.df[res_ctx_nTg_vs_pSyn_pos.df$padj < 0.05, ]
nTg_vs_pSyn_pos_SIGS.df <- na.omit(nTg_vs_pSyn_pos_SIGS.df)

nTg_vs_pSyn_neg_SIGS.df <- res_ctx_nTg_vs_pSyn_neg.df[res_ctx_nTg_vs_pSyn_neg.df$padj < 0.05, ]
nTg_vs_pSyn_neg_SIGS.df <- na.omit(nTg_vs_pSyn_neg_SIGS.df)


gene_counts1 <- table(nTg_vs_pSyn_pos_SIGS.df$gene)
genes_keep1 <- names(gene_counts1[gene_counts1 >= 2])

gene_counts2 <- table(nTg_vs_pSyn_neg_SIGS.df$gene)
genes_keep2 <- names(gene_counts2[gene_counts2 >= 2])

genes_keep <- union(genes_keep1, genes_keep2)

genes_for_viz <- as.data.frame(genes_keep)

genes_to_keep <- read.csv2("/xenium_Line61_final/xenium_upload1/general/genes_to_keep_DE.csv", header = F)
genes_to_keep <- genes_to_keep$V1

genes_for_viz <- genes_for_viz[genes_for_viz$genes_keep %in% genes_to_keep, ]
saveRDS(genes_for_viz, file = "/xenium_Line61_final/ctx/outs/DE_genes_for_viz.rds")


################################################################################

# filtering results for genes, making matrices for DE (log2fc / FDR)

res_ctx_nTg_vs_pSyn_pos.df <- readRDS("/xenium_Line61_final/ctx/outs/DE_ctx_nTg_vs_pSyn_pos_ALL.rds")
res_ctx_nTg_vs_pSyn_neg.df <- readRDS("/xenium_Line61_final/ctx/outs/DE_ctx_nTg_vs_pSyn_neg_ALL.rds")
res_Tg_pos_vs_neg.df <- readRDS("/xenium_Line61_final/ctx/outs/DE_ctx_Tg_only_ALL.rds")
genes_for_viz <- readRDS("/xenium_Line61_final/ctx/outs/DE_genes_for_viz.rds")

nTg_vs_pSyn_pos_VIZ.df <- res_ctx_nTg_vs_pSyn_pos.df[res_ctx_nTg_vs_pSyn_pos.df$gene %in% genes_for_viz, ]
nTg_vs_pSyn_neg_VIZ.df <- res_ctx_nTg_vs_pSyn_neg.df[res_ctx_nTg_vs_pSyn_neg.df$gene %in% genes_for_viz, ]
Tg_pos_vs_neg_VIZ.df <- res_Tg_pos_vs_neg.df[res_Tg_pos_vs_neg.df$gene %in% genes_for_viz, ]


###
l2fc_pSyn_pos <- nTg_vs_pSyn_pos_VIZ.df[c("gene", "celltype", "log2FoldChange")]
l2fc_pSyn_pos <- pivot_wider(l2fc_pSyn_pos, names_from = celltype, values_from = log2FoldChange)
saveRDS(l2fc_pSyn_pos, file = "/xenium_Line61_final/ctx/outs/l2fc_nTg_Tg_pSyn_pos_VIZ.rds")

fdr_pSyn_pos <- nTg_vs_pSyn_pos_VIZ.df[c("gene", "celltype", "padj")]
fdr_pSyn_pos <- pivot_wider(fdr_pSyn_pos, names_from = celltype, values_from = padj)
saveRDS(fdr_pSyn_pos, file = "/xenium_Line61_final/ctx/outs/fdr_nTg_Tg_pSyn_pos_VIZ.rds")


###
l2fc_pSyn_neg <- nTg_vs_pSyn_neg_VIZ.df[c("gene", "celltype", "log2FoldChange")]
l2fc_pSyn_neg <- pivot_wider(l2fc_pSyn_neg, names_from = celltype, values_from = log2FoldChange)
saveRDS(l2fc_pSyn_neg, file = "/xenium_Line61_final/ctx/outs/l2fc_nTg_Tg_pSyn_neg_VIZ.rds")

fdr_pSyn_neg <- nTg_vs_pSyn_neg_VIZ.df[c("gene", "celltype", "padj")]
fdr_pSyn_neg <- pivot_wider(fdr_pSyn_neg, names_from = celltype, values_from = padj)
saveRDS(fdr_pSyn_neg, file = "/xenium_Line61_final/ctx/outs/fdr_nTg_Tg_pSyn_neg_VIZ.rds")


###
l2fc_Tg_only <- Tg_pos_vs_neg_VIZ.df[c("gene", "celltype", "log2FoldChange")]
l2fc_Tg_only <- pivot_wider(l2fc_Tg_only, names_from = celltype, values_from = log2FoldChange)
saveRDS(l2fc_Tg_only, file = "/xenium_Line61_final/ctx/outs/l2fc_Tg_only_VIZ.rds")

fdr_Tg_only <- Tg_pos_vs_neg_VIZ.df[c("gene", "celltype", "padj")]
fdr_Tg_only <- pivot_wider(fdr_Tg_only, names_from = celltype, values_from = padj)
saveRDS(fdr_Tg_only, file = "/xenium_Line61_final/ctx/outs/fdr_Tg_only_VIZ.rds")

#######

nTg_vs_pSyn_pos_CUSTOM.df <- res_ctx_nTg_vs_pSyn_pos.df[res_ctx_nTg_vs_pSyn_pos.df$gene %in% genes_to_keep, ]
nTg_vs_pSyn_neg_CUSTOM.df <- res_ctx_nTg_vs_pSyn_neg.df[res_ctx_nTg_vs_pSyn_neg.df$gene %in% genes_to_keep, ]
Tg_pos_vs_neg_CUSTOM.df <- res_Tg_pos_vs_neg.df[res_Tg_pos_vs_neg.df$gene %in% genes_to_keep, ]


write.csv(nTg_vs_pSyn_pos_CUSTOM.df, file = "/xenium_Line61_final/ctx/outs/DE_nTg_vs_pSyn_pos_CUSTOM.csv")
write.csv(nTg_vs_pSyn_neg_CUSTOM.df, file = "/xenium_Line61_final/ctx/outs/DE_nTg_vs_pSyn_neg_CUSTOM.csv")
write.csv(Tg_pos_vs_neg_CUSTOM.df, file = "/xenium_Line61_final/ctx/outs/DE_pSyn_pos_vs_neg_CUSTOM.csv")
