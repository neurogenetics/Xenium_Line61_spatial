library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(extrafont)
library(svglite)
library(monocle3)
library(SeuratWrappers)


xenium_ctx_tg <- readRDS("/xenium_Line61_final/ctx/xenium_CTX_tg_only_for_DESeq.rds")

genes_for_viz <- readRDS("/xenium_Line61_final/ctx/outs/DE_genes_for_viz.rds")


################################################################################

# subset for just cortical ExNs

Idents(xenium_ctx_tg) <- 'celltype'
ExNs <- WhichCells(xenium_ctx_tg, idents = c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", "L5 ET ExN", 
                                             "L5 NP ExN", "L6 IT ExN", "L6 CT ExN", "L6b ExN")) 

ctx_ExNs <- subset(xenium_ctx_tg, cells = ExNs)


# calculate hSNCA expression and assign in metadata

orig_ident <- ctx_ExNs@meta.data$orig.ident
has_pSyn <- ctx_ExNs@meta.data$has_pSyn
celltype <- ctx_ExNs@meta.data$celltype

expression_matrix <- as.data.frame(LayerData(ctx_ExNs, assay = "Xenium", layer = "data"))
hSNCA_expression <- expression_matrix["hSNCA", ]
hSNCA_expression <- t(hSNCA_expression)

df <- data.frame(hSNCA_expression = hSNCA_expression)

df$hSNCA <- as.numeric(df$hSNCA)
df$expm1 <- expm1(df$hSNCA)

expm1 <- df$expm1

ctx_ExNs <- AddMetaData(ctx_ExNs, metadata = expm1, col.name = "hSNCA_exp")


# add bins for hSNCA expression

ctx_ExNs$hSNCA_bin <- ifelse(ctx_ExNs$hSNCA_exp <= 100, "0-100",
                             ifelse(ctx_ExNs$hSNCA_exp > 100 & ctx_ExNs$hSNCA_exp <= 200, "100-200",
                                    ifelse(ctx_ExNs$hSNCA_exp > 200 & ctx_ExNs$hSNCA_exp <= 300, "200-300",
                                           ifelse(ctx_ExNs$hSNCA_exp > 300 & ctx_ExNs$hSNCA_exp <= 400, "300-400",
                                                  ifelse(ctx_ExNs$hSNCA_exp > 400 & ctx_ExNs$hSNCA_exp <= 500, "400-500",
                                                         ifelse(ctx_ExNs$hSNCA_exp > 500 & ctx_ExNs$hSNCA_exp <= 600, "500-600",
                                                                ifelse(ctx_ExNs$hSNCA_exp > 600 & ctx_ExNs$hSNCA_exp <= 700, "600-700",
                                                                       ifelse(ctx_ExNs$hSNCA_exp > 700 & ctx_ExNs$hSNCA_exp <= 800, "700-800",
                                                                              ifelse(ctx_ExNs$hSNCA_exp > 800 & ctx_ExNs$hSNCA_exp <= 900, "800-900",
                                                                                     ifelse(ctx_ExNs$hSNCA_exp > 900 & ctx_ExNs$hSNCA_exp <= 1000, "900-1000", "1000+"))))))))))

ctx_ExNs$hSNCA_bin <- factor(ctx_ExNs$hSNCA_bin, levels = c("0-100", "100-200", "200-300", "300-400", "400-500", 
                                                            "500-600", "600-700", "700-800", "800-900", "900-1000", "1000+"))

################################################################################

L23_IT <- subset(ctx_ExNs, cells = WhichCells(ctx_ExNs, idents = "L2/3 IT ExN"))
L45_IT <- subset(ctx_ExNs, cells = WhichCells(ctx_ExNs, idents = "L4/5 IT ExN"))
L5_IT <- subset(ctx_ExNs, cells = WhichCells(ctx_ExNs, idents = "L5 IT ExN"))
L5_ET <- subset(ctx_ExNs, cells = WhichCells(ctx_ExNs, idents = "L5 ET ExN"))
L5_NP <- subset(ctx_ExNs, cells = WhichCells(ctx_ExNs, idents = "L5 NP ExN"))
L6_IT <- subset(ctx_ExNs, cells = WhichCells(ctx_ExNs, idents = "L6 IT ExN"))
L6_CT <- subset(ctx_ExNs, cells = WhichCells(ctx_ExNs, idents = "L6 CT ExN"))
L6b <- subset(ctx_ExNs, cells = WhichCells(ctx_ExNs, idents = "L6b ExN"))

cds_list <- c(L23_IT, L45_IT, L5_IT, L5_ET, L5_NP, L6_IT, L6_CT, L6b)

hSNCA_terms_list <- list()

for (i in seq_along(cds_list)) {
  cds <- cds_list[[i]]
  cds <- SeuratWrappers::as.cell_data_set(cds)
  cds <- preprocess_cds(cds, num_dim = 30)
  cds <- reduce_dimension(cds)
  gene_fits <- fit_models(cds, model_formula_str = "~hSNCA_exp")
  coefs <- coefficient_table(gene_fits)
  hSNCA_terms <- coefs %>% filter(term == "hSNCA_exp")
  hSNCA_terms <- hSNCA_terms %>%
    filter(gene_id %in% genes_for_viz) %>%
    select(gene_id, q_value, estimate)
  
  hSNCA_terms_list[[i]] <- hSNCA_terms
}


# generate dfs with the q-values and estimates for complexheatmap

estimate_list <- list()
q_value_list <- list()

for (i in seq_along(hSNCA_terms_list)) {
  hSNCA_terms <- hSNCA_terms_list[[i]]
  
  estimate_list[[i]] <- hSNCA_terms %>% select(gene_id, estimate)
  q_value_list[[i]] <- hSNCA_terms %>% select(gene_id, q_value)
}

# combine the lists into dataframes
estimate.df <- do.call(cbind, lapply(estimate_list, function(x) setNames(x$estimate, x$gene_id)))
q_value.df <- do.call(cbind, lapply(q_value_list, function(x) setNames(x$q_value, x$gene_id)))

estimate.df <- as.data.frame(estimate.df)
q_value.df <- as.data.frame(q_value.df)

rownames(estimate.df) <- estimate_list[[1]]$gene_id
rownames(q_value.df) <- q_value_list[[1]]$gene_id


colnames(estimate.df) <- c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", "L5 ET ExN", 
                           "L5 NP ExN", "L6 IT ExN", "L6 CT ExN", "L6b ExN")

colnames(q_value.df) <- c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", "L5 ET ExN", 
                          "L5 NP ExN", "L6 IT ExN", "L6 CT ExN", "L6b ExN")


# multiply values in estimate.df by 10^4
estimate.df <- estimate.df * 10^4


saveRDS(estimate.df, file = "/xenium_Line61_final/ctx/outs/GLM_estimates.rds")
saveRDS(q_value.df, file = "/xenium_Line61_final/ctx/outs/GLM_qvals.rds")

write.csv(estimate.df, "/xenium_Line61_final/ctx/outs/GLM_estimates.csv")
write.csv(q_value.df, "/xenium_Line61_final/ctx/outs/GLM_qvals.csv")
