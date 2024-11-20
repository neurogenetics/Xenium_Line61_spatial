library(tidyverse)
library(dplyr)
library(eulerr)
library(svglite)


# filtering results for genes, making matrices for DE (log2fc / FDR)

res_ctx_nTg_vs_pSyn_pos.df <- readRDS("/xenium_Line61_final/ctx/outs/DE_ctx_nTg_vs_pSyn_pos_ALL.rds")
res_ctx_nTg_vs_pSyn_neg.df <- readRDS("/xenium_Line61_final/ctx/outs/DE_ctx_nTg_vs_pSyn_neg_ALL.rds")
res_Tg_pos_vs_neg.df <- readRDS("/xenium_Line61_final/ctx/outs/DE_ctx_Tg_only_ALL.rds")

xenium_genes.df <- as.data.frame(read.csv2("/xenium_Line61_final/xenium_upload1/general/xenium_genes.csv", 
                                           ",", header = T))
all_genes <- xenium_genes.df$gene_all
custom_genes <- all_genes[249:347]


# list of dataframes
dfs <- list(res_ctx_nTg_vs_pSyn_pos.df, res_ctx_nTg_vs_pSyn_neg.df, res_Tg_pos_vs_neg.df)

# filter significant DEGs and split by cell type
filtered_dfs <- lapply(dfs, function(df) {
  df <- df[df$gene %in% custom_genes, ]
  df <- df[df$padj <= 0.05, ]
  split(df, df$celltype)
})

# get unique cell types 
cell_types <- unique(res_ctx_nTg_vs_pSyn_pos.df$celltype)

# create the final list for each cell type
final_list <- lapply(cell_types, function(cell_type) {
  lapply(filtered_dfs, function(filtered_df) {
    if (cell_type %in% names(filtered_df)) {
      filtered_df[[cell_type]]$gene
    } 
  })
})

# name the final list with cell types
names(final_list) <- cell_types


# split list into cell types
list2env(final_list, envir = .GlobalEnv)


names(`L2/3 IT ExN`) <- c("non-tg vs. pSyn+ α-syn-tg", 
                          "non-tg vs. pSyn- α-syn-tg", 
                          "pSyn+ vs. pSyn- α-syn-tg")

names(`L4/5 IT ExN`) <- c("non-tg vs. pSyn+ α-syn-tg", 
                          "non-tg vs. pSyn- α-syn-tg", 
                          "pSyn+ vs. pSyn- α-syn-tg")

names(`L5 IT ExN`) <- c("non-tg vs. pSyn+ α-syn-tg", 
                        "non-tg vs. pSyn- α-syn-tg", 
                        "pSyn+ vs. pSyn- α-syn-tg")

names(`L5 ET ExN`) <- c("non-tg vs. pSyn+ α-syn-tg", 
                        "non-tg vs. pSyn- α-syn-tg", 
                        "pSyn+ vs. pSyn- α-syn-tg")

names(`L5 NP ExN`) <- c("non-tg vs. pSyn+ α-syn-tg", 
                        "non-tg vs. pSyn- α-syn-tg", 
                        "pSyn+ vs. pSyn- α-syn-tg")

names(`L6 IT ExN`) <- c("non-tg vs. pSyn+ α-syn-tg", 
                        "non-tg vs. pSyn- α-syn-tg", 
                        "pSyn+ vs. pSyn- α-syn-tg")

names(`L6 CT ExN`) <- c("non-tg vs. pSyn+ α-syn-tg", 
                        "non-tg vs. pSyn- α-syn-tg", 
                        "pSyn+ vs. pSyn- α-syn-tg")

names(`L6b ExN`) <- c("non-tg vs. pSyn+ α-syn-tg", 
                      "non-tg vs. pSyn- α-syn-tg", 
                      "pSyn+ vs. pSyn- α-syn-tg")

################################################################################

# making euler plots (Supplementary Fig. 7a)

fit1 <- euler(`L2/3 IT ExN`)
pdf("/xenium_Line61_final/ctx/outs/DE_euler_L23_IT.pdf", 
    width = 2.5, height = 2.5, family = "Helvetica")
plot(fit1, 
     quantities = TRUE,
     fills = list(fill = c("#c26a77", "steelblue4", "#dcce7d"), alpha = 0.8),
     labels = NULL) 
dev.off()



fit2 <- euler(`L4/5 IT ExN`)
pdf("/xenium_Line61_final/ctx/outs/DE_euler_L45_IT.pdf", 
    width = 2.5, height = 2.5, family = "Helvetica")
plot(fit2, 
     quantities = TRUE,
     fills = list(fill = c("#c26a77", "steelblue4", "#dcce7d"), alpha = 0.8),
     labels = NULL) 
dev.off()



fit3 <- euler(`L5 IT ExN`)
pdf("/xenium_Line61_final/ctx/outs/DE_euler_L5_IT.pdf", 
    width = 2.5, height = 2.5, family = "Helvetica")
plot(fit3, 
     quantities = TRUE,
     fills = list(fill = c("#c26a77", "steelblue4", "#dcce7d"), alpha = 0.8),
     labels = NULL) 
dev.off()



fit4 <- euler(`L5 ET ExN`)
pdf("/xenium_Line61_final/ctx/outs/DE_euler_L5_ET.pdf", 
    width = 2.5, height = 2.5, family = "Helvetica")
plot(fit4, 
     quantities = TRUE,
     fills = list(fill = c("#c26a77", "steelblue4", "#dcce7d"), alpha = 0.8),
     labels = NULL) 
dev.off()



fit5 <- euler(`L5 NP ExN`)
pdf("/xenium_Line61_final/ctx/outs/DE_euler_L5_NP.pdf", 
    width = 2.5, height = 2.5, family = "Helvetica")
plot(fit5, 
     quantities = TRUE,
     fills = list(fill = c("#c26a77", "steelblue4", "#dcce7d"), alpha = 0.8),
     labels = NULL) 
dev.off()



fit6 <- euler(`L6 IT ExN`)
pdf("/xenium_Line61_final/ctx/outs/DE_euler_L6_IT.pdf", 
    width = 2.5, height = 2.5, family = "Helvetica")
plot(fit6, 
     quantities = TRUE,
     fills = list(fill = c("#c26a77", "steelblue4", "#dcce7d"), alpha = 0.8),
     labels = NULL) 
dev.off()



fit7 <- euler(`L6 CT ExN`)
pdf("/xenium_Line61_final/ctx/outs/DE_euler_L6_CT.pdf", 
    width = 2.5, height = 2.5, family = "Helvetica")
plot(fit7, 
     quantities = TRUE,
     fills = list(fill = c("#c26a77", "steelblue4", "#dcce7d"), alpha = 0.8),
     labels = NULL) 
dev.off()



fit8 <- euler(`L6b ExN`)
pdf("/xenium_Line61_final/ctx/outs/DE_euler_L6b.pdf", 
    width = 2.5, height = 2.5, family = "Helvetica")
plot(fit8, 
     quantities = TRUE,
     fills = list(fill = c("#c26a77", "steelblue4", "#dcce7d"), alpha = 0.8),
     labels = NULL) 
dev.off()
