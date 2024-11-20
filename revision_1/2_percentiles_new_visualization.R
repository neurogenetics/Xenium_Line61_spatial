library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(broom)
library(ggrepel)
library(extrafont)
library(svglite)

xenium_ctx_tg <- readRDS("/xenium_Line61_final/ctx/xenium_CTX_tg_only_for_DESeq.rds")

################################################################################

# generate dataframe w/ hSNCA expression and pSyn status + percentiles

orig_ident <- xenium_ctx_tg@meta.data$orig.ident
has_pSyn <- xenium_ctx_tg@meta.data$has_pSyn
celltype <- xenium_ctx_tg@meta.data$celltype

expression_matrix <- as.data.frame(LayerData(xenium_ctx_tg, assay = "Xenium", layer = "data"))
hSNCA_expression <- expression_matrix["hSNCA", ]
hSNCA_expression <- t(hSNCA_expression)

df <- data.frame(Cell_ID = rownames(hSNCA_expression),
                 orig.ident = orig_ident,
                 has.pSyn = has_pSyn,
                 celltype = celltype,
                 hSNCA_expression = hSNCA_expression)

celltypes_keep <- c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", "L5 ET ExN", "L5 NP ExN", 
                    "L6 IT ExN", "L6 CT ExN", "L6b ExN")

df <- df[df$celltype %in% celltypes_keep, ]
df$hSNCA <- as.numeric(df$hSNCA)
df$expm1 <- expm1(df$hSNCA)

df$celltype <- factor(df$celltype, levels = c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", "L5 ET ExN", "L5 NP ExN", 
                                              "L6 IT ExN", "L6 CT ExN", "L6b ExN"))

pSyn_pos.df <- as.data.frame(df[df$has.pSyn %in% "TRUE", ])

percentiles1 <- aggregate(expm1 ~ orig.ident + celltype, data = pSyn_pos.df, FUN = function(x) {
  c(p25 = quantile(x, probs = 0.25))
})

# top 25th percentile for pSyn-

pSyn_neg.df <- df[df$has.pSyn %in% "FALSE", ]

percentiles2 <- aggregate(expm1 ~ orig.ident + celltype, data = pSyn_neg.df, FUN = function(x) {
  c(p75 = quantile(x, probs = 0.75))
})

merged.df <- df %>%
  left_join(percentiles1, by = c("orig.ident", "celltype")) %>%
  left_join(percentiles2, by = c("orig.ident", "celltype"))

colnames(merged.df) <- c("Cell_ID", "orig.ident", "has.pSyn", "celltype", "hSNCA", "expm1", "p25_pos", "p75_neg")

merged.df$expm1 <- as.numeric(merged.df$expm1)
merged.df$p25_pos <- as.numeric(merged.df$p25_pos)
merged.df$p75_neg <- as.numeric(merged.df$p75_neg)

merged.df$group <- ifelse(merged.df$has.pSyn == TRUE & merged.df$expm1 <= merged.df$p25_pos, "pSyn+/hSNCA-low",
                          ifelse(merged.df$has.pSyn == FALSE & merged.df$expm1 >= merged.df$p75_neg, "pSyn-/hSNCA-high",
                                 ifelse(merged.df$has.pSyn == FALSE & merged.df$expm1 <= merged.df$p75_neg, "pSyn-/hSNCA-low", "pSyn+/hSNCA-high")))
merged.df$group <- factor(merged.df$group, 
                          levels = c("pSyn-/hSNCA-low", "pSyn-/hSNCA-high", "pSyn+/hSNCA-low", "pSyn+/hSNCA-high"))


p1 <- ggplot(merged.df, aes(x = celltype, y = expm1, color = has.pSyn)) +
  geom_boxplot(outliers = F, coef = 0) +
  geom_point(position = position_dodge(width = 0.75), aes(fill = group), shape = 21, size = 3) +
  scale_color_discrete(name = "pSyn status",
                       breaks = c(TRUE, FALSE)) +
  scale_color_manual(values = c("black", "black"),
                     labels = c("pSyn-", "pSyn+")) + 
  scale_fill_manual(values = c("pSyn+/hSNCA-low" = "#fab2ac", "pSyn-/hSNCA-high" = "#035c80",
                               "pSyn+/hSNCA-high" = "darkred", "pSyn-/hSNCA-low" = "#f0f7fa")) +
  labs(x = "",
       y = "hSNCA expression") + 
  theme_bw() +
  theme(text = element_text(family = "Arial"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold"),
        legend.title = element_blank(), 
        legend.text = element_text(size = 10, face = "bold")) 

svglite("/xenium_Line61_final/ctx/outs/REVISED_ctx_hSNCA_percentiles.svg", 
        width = 14.5, height = 4.3)
plot(p1)
dev.off()
