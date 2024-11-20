library(Seurat)
library(DESeq2)
library(tidyverse)
library(dplyr)
library(magrittr)

xenium_hip <- readRDS("/xenium_Line61_final/hip/xenium_HIP_metadata_pSyn.rds")

################################################################################

# retaining Tg cells (for pSyn+ vs. pSyn-)

Idents(xenium_hip) <- "genotype"
tg_cells_hip <- WhichCells(xenium_hip, idents = "Tg")
# using subset_opt function found here: https://github.com/alikhuseynov/add-on_R/blob/develop/R/subset_obj_seurat_v2.R 
# because Seurat is bugged and won't let you subset if there aren't some cells from every fov in the subset object
xenium_hip_tg <- subset_opt(xenium_hip, cells = tg_cells_hip) 
saveRDS(xenium_hip_tg, file = "/xenium_Line61_final/hip/xenium_HIP_tg_only_for_DESeq.rds")

################################################################################

hip_ExNs <- c("CA1", "CA2/3", "DG")

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

# Tg pSyn+ vs. pSyn-
  
xenium_tg_only <- readRDS("/xenium_Line61_final/hip/xenium_HIP_tg_only_for_DESeq.rds")

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


res_Tg_pos_vs_neg <- run_DESeq_Tg_only(cts3_split_modified, hip_ExNs)
res_Tg_pos_vs_neg$genes <- rownames(res_Tg_pos_vs_neg)
res_Tg_pos_vs_neg$genes <- gsub("^result_", "", res_Tg_pos_vs_neg$genes)
res_Tg_pos_vs_neg <- separate(res_Tg_pos_vs_neg, genes, into = c("celltype", "gene"), sep = "\\.")

res_Tg_pos_vs_neg.df <- as.data.frame(res_Tg_pos_vs_neg)

saveRDS(res_Tg_pos_vs_neg.df, file = "/xenium_Line61_final/hip/outs/DE_hip_Tg_only_ALL.rds")

################################################################################

genes_to_keep <- read.csv2("/xenium_Line61_final/xenium_upload1/general/genes_to_keep_DE.csv", header = F)
genes_to_keep <- genes_to_keep$V1

hip_DE <- res_Tg_pos_vs_neg.df[res_Tg_pos_vs_neg.df$gene %in% genes_to_keep, ]

hip_DE <- hip_DE[hip_DE$padj < 0.05, ]
hip_DE <- na.omit(hip_DE)
