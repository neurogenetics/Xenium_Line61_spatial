library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(broom)
library(ggrepel)
library(extrafont)
library(svglite)

pSyn_counts.df <- readRDS("/xenium_Line61_final/ctx/outs/pSyn_counts.rds")

AvgExp_ctx.df <- readRDS("/xenium_Line61_final/ctx/outs/avgexp_ctx_df.rds")

################################################################################


pSyn_avg <- pSyn_counts.df %>%
  group_by(cell_types) %>%
  summarize(mean = mean(pct))

pSyn_avg.df <- as.data.frame(pSyn_avg)


################################################################################
################################################################################
################################################################################

# for analyzing non-tg samples only

nTg_keep <- c("nTg.1", "nTg.2", "nTg.3", "nTg.4")

AvgExp_nTg.df <- AvgExp_ctx.df[AvgExp_ctx.df$animal %in% nTg_keep, ]

AvgExp_nTg.df <- as.data.frame(AvgExp_nTg.df %>%
                                 group_by(gene, celltype) %>%
                                 summarize(mean = mean(Expression)))

nTg_merged.df <- merge(AvgExp_nTg.df, pSyn_avg.df, by.x = "celltype", by.y = "cell_types")

colnames(nTg_merged.df) <- c("celltype", "gene_name", "expression", "pSyn_avg")

################################################################################

# function for iteratively running linear regressions b/w gene expression and % pSyn+ 

calculate_gene_correlation <- function(df) {
  gene_list <- unique(df$gene_name)
  
  results <- lapply(gene_list, function(gene) {
    subset_df <- df %>% filter(gene_name == gene)
    
    lm_result <- lm(pSyn_avg ~ expression, data = subset_df)
    tidy_result <- tidy(lm_result)
    
    r_squared <- summary(lm_result)$r.squared
    p_value <- tidy_result$`p.value`[2]
    
    return(data.frame(gene = gene, r_squared = r_squared, p_value = p_value))
  })
  
  results_df <- do.call(rbind, results)
  return(results_df)
}

################################################################################

# run fxn for non-tg samples (all cell types)

nTg_linreg.df <- as.data.frame(calculate_gene_correlation(nTg_merged.df))

################################################################################

# graphing Plk2 correlation (all cell types) (Fig. 2f)

plk2_all.df <- nTg_merged.df[nTg_merged.df$gene_name == "Plk2", ]
rownames(plk2_all.df) <- NULL
plk2_all.df <- column_to_rownames(plk2_all.df, "celltype")

p1 <- ggplot(plk2_all.df, aes(x = expression, y = pSyn_avg)) + 
  geom_point(size = 2) + 
  geom_smooth(method='lm') + 
  geom_label_repel(aes(expression, pSyn_avg, label = rownames(plk2_all.df))) + 
  labs(x = "Mean Plk2 expression, non-tg only (expm1[log-norm counts])",
       y = "% cells pSyn+") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        text = element_text(family = "Arial")) +
  annotate("text", x = 20, y = 56, label = "r^2 = 0.7051") + 
  annotate("text", x = 20, y = 53, label = "p = 9.00e-05")

svglite("/xenium_Line61_final/ctx/outs/ctx_plk2_linreg.svg", 
        width = 9.5, height = 3.8)
plot(p1)
dev.off()

################################################################################

# analyzing correlations w/ only ExNs

ctx_ExNs <- c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", "L5 ET ExN", 
              "L5 NP ExN", "L6 IT ExN", "L6 CT ExN", "L6b ExN")

nTg_ExNs_merged.df <- nTg_merged.df[nTg_merged.df$celltype %in% ctx_ExNs, ]

nTg_linreg_ExNs.df <- as.data.frame(calculate_gene_correlation(nTg_ExNs_merged.df))


################################################################################
################################################################################
################################################################################

# for analyzing Tg samples only

Tg_keep <- c("Tg.1", "Tg.2", "Tg.3", "Tg.4")

AvgExp_Tg.df <- AvgExp_ctx.df[AvgExp_ctx.df$animal %in% Tg_keep, ]

AvgExp_Tg.df <- as.data.frame(AvgExp_Tg.df %>%
                                group_by(gene, celltype) %>%
                                summarize(mean = mean(Expression)))

Tg_merged.df <- merge(AvgExp_Tg.df, pSyn_avg.df, by.x = "celltype", by.y = "cell_types")

colnames(Tg_merged.df) <- c("celltype", "gene_name", "expression", "pSyn_avg")

################################################################################

# run fxn for all Tg cell types

Tg_linreg.df <- as.data.frame(calculate_gene_correlation(Tg_merged.df))

################################################################################

# graphing hSNCA correlation (all cell types)

hsnca_all.df <- Tg_merged.df[Tg_merged.df$gene_name == "hSNCA", ]
rownames(hsnca_all.df) <- NULL
hsnca_all.df <- column_to_rownames(hsnca_all.df, "celltype")

p2 <- ggplot(hsnca_all.df, aes(x = expression, y = pSyn_avg)) + 
  geom_point(size = 2) + 
  geom_smooth(method='lm') + 
  geom_label_repel(aes(expression, pSyn_avg, label = rownames(hsnca_all.df))) + 
  labs(x = "Mean hSNCA expression (expm1[log-norm counts])",
       y = "% cells pSyn+") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        text = element_text(family = "Arial")) +
  annotate("text", x = 250, y = 56, label = "r^2 = 0.7702") + 
  annotate("text", x = 250, y = 53, label = "p = 1.71e-05")

svglite("/xenium_Line61_final/ctx/outs/ctx_hSNCA_linreg_allcelltypes.svg", 
        width = 9.5, height = 3.8)
plot(p2)
dev.off()

################################################################################

# analyzing only Tg ExNs

ctx_ExNs <- c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", "L5 ET ExN", 
              "L5 NP ExN", "L6 IT ExN", "L6 CT ExN", "L6b ExN")

Tg_ExNs_merged.df <- Tg_merged.df[Tg_merged.df$celltype %in% ctx_ExNs, ]

Tg_linreg_ExNs.df <- as.data.frame(calculate_gene_correlation(Tg_ExNs_merged.df))

################################################################################

hsnca_ExNs.df <- Tg_ExNs_merged.df[Tg_ExNs_merged.df$gene_name == "hSNCA", ]
rownames(hsnca_ExNs.df) <- NULL
hsnca_ExNs.df <- column_to_rownames(hsnca_ExNs.df, "celltype")

p3 <- ggplot(hsnca_ExNs.df, aes(x = expression, y = pSyn_avg)) + 
  geom_point(size = 2) + 
  geom_smooth(method='lm') + 
  geom_label_repel(aes(expression, pSyn_avg, label = rownames(hsnca_ExNs.df))) + 
  labs(x = "Mean hSNCA expression (expm1[log-norm counts])",
       y = "% cells pSyn+") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        text = element_text(family = "Arial")) +
  annotate("text", x = 375, y = 56, label = "r^2 = 0.3635") + 
  annotate("text", x = 375, y = 53, label = "p = 0.1136")

svglite("/xenium_Line61_final/ctx/outs/ctx_hSNCA_linreg_ExNs.svg", 
        width = 9.5, height = 3.8)
plot(p3)
dev.off()

################################################################################

write.csv(nTg_linreg.df, file = "/xenium_Line61_final/ctx/outs/nTg_linreg.csv")

write.csv(nTg_linreg_ExNs.df, file = "/xenium_Line61_final/ctx/outs/nTg_linreg_ExNs.csv")

write.csv(Tg_linreg.df, file = "/xenium_Line61_final/ctx/outs/Tg_linreg.csv")

write.csv(Tg_linreg_ExNs.df, file = "/xenium_Line61_final/ctx/outs/Tg_linreg_ExNs.csv")
