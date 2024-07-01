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
xenium_ctx <- readRDS("/xenium_Line61_final/ctx/xenium_CTX_metadata_pSyn.rds")

################################################################################

# plot showing individual Tg cells with their hSNCA expression, split by pSyn+/- (Supplementary Fig. 6a)

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

p1 <- ggplot(df, aes(x = celltype, y = expm1, color = has.pSyn)) +
  geom_boxplot(outliers = F, coef = 0) +
  geom_point(position = position_dodge(width = 0.75)) +
  scale_color_discrete(name = "pSyn status",
                       breaks = c(TRUE, FALSE)) +
  scale_color_manual(values = c("grey", "red"),
                     labels = c("pSyn-", "pSyn+")) +
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

svglite("/xenium_Line61_final/ctx/outs/ctx_hSNCA_exp_by_celltype.svg", 
        width = 10.5, height = 4)
plot(p1)
dev.off()

################################################################################
################################################################################
################################################################################

# extracting cells in the top 25% by hSNCA exp for pSyn- cells, and bottom 25% by hSNCA exp for pSyn+ cells

################################################################################

# bottom 25th percentile for pSyn+

pSyn_pos.df <- as.data.frame(df[df$has.pSyn %in% "TRUE", ])

percentiles1 <- aggregate(expm1 ~ orig.ident + celltype, data = pSyn_pos.df, FUN = function(x) {
  c(p25 = quantile(x, probs = 0.25))
})

# top 25th percentile for pSyn-

pSyn_neg.df <- df[df$has.pSyn %in% "FALSE", ]

percentiles2 <- aggregate(expm1 ~ orig.ident + celltype, data = pSyn_neg.df, FUN = function(x) {
  c(p75 = quantile(x, probs = 0.75))
})

# combine dfs

merged.df <- df %>%
  left_join(percentiles1, by = c("orig.ident", "celltype")) %>%
  left_join(percentiles2, by = c("orig.ident", "celltype"))

colnames(merged.df) <- c("Cell_ID", "orig.ident", "has.pSyn", "celltype", "hSNCA", "expm1", "p25_pos", "p75_neg")

merged.df$expm1 <- as.numeric(merged.df$expm1)
merged.df$p25_pos <- as.numeric(merged.df$p25_pos)
merged.df$p75_neg <- as.numeric(merged.df$p75_neg)


################################################################################


# extract cell IDs of cells within the defined parameters (top 25% pSyn neg or bottom 25% pSyn pos)

selected_cells <- merged.df %>%
  filter(
    (has.pSyn == TRUE & expm1 <= p25_pos) |  
      (has.pSyn == FALSE & expm1 >= p75_neg)   
  ) %>%
  select(Cell_ID)  

cells_keep <- selected_cells$Cell_ID


# subset object for just tg cells in the percentiles

xenium_ctx_percentiles <- subset(xenium_ctx_tg, cells = cells_keep)
saveRDS(xenium_ctx_percentiles, file = "/xenium_Line61_final/ctx/xenium_CTX_tg_percentiles.rds")



# subset object for non-tg cells and tg cells in the percentiles

Idents(xenium_ctx) <- "genotype"
ntg_keep <- WhichCells(xenium_ctx, idents = "nTg")

all_cells_keep <- c(ntg_keep, cells_keep)

xenium_ctx_percentiles_w_nTg <- subset(xenium_ctx, cells = all_cells_keep)
saveRDS(xenium_ctx_percentiles_w_nTg, 
        file = "/xenium_Line61_final/ctx/xenium_CTX_tg_percentiles_w_nTg.rds")


################################################################################

# graph to show 25th and 75th percentiles (Fig. 5b)

xenium_ctx_percentiles <- readRDS("/xenium_Line61_final/ctx/xenium_CTX_tg_percentiles.rds")

orig_ident <- xenium_ctx_percentiles@meta.data$orig.ident
has_pSyn <- xenium_ctx_percentiles@meta.data$has_pSyn
celltype <- xenium_ctx_percentiles@meta.data$celltype

expression_matrix <- as.data.frame(LayerData(xenium_ctx_percentiles, assay = "Xenium", layer = "data"))
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

p2 <- ggplot(df, aes(x = celltype, y = expm1, color = has.pSyn)) +
  geom_boxplot(outliers = F, coef = 0) +
  geom_point(position = position_dodge(width = 0.75)) +
  scale_color_discrete(name = "pSyn status",
                       breaks = c(TRUE, FALSE)) +
  scale_color_manual(values = c("grey", "red"),
                     labels = c("pSyn-", "pSyn+")) +
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

svglite("/xenium_Line61_final/ctx/outs/ctx_hSNCA_exp_by_celltype_percentiles.svg", 
        width = 10.5, height = 4)
plot(p2)
dev.off()

################################################################################

# getting cell IDs from different percentiles for visualization in Xenium Explorer (for Fig. 5a)


# percentiles for pSyn+ cells

pSyn_pos.df <- as.data.frame(df[df$has.pSyn %in% "TRUE", ])

percentiles_pos <- aggregate(expm1 ~ orig.ident + celltype, data = pSyn_pos.df, FUN = function(x) {
  c(p25 = quantile(x, probs = 0.25),
    p50 = quantile(x, probs = 0.50),
    p75 = quantile(x, probs = 0.75))
})

merged_pos.df <- merge(pSyn_pos.df, percentiles_pos, by = c("orig.ident", "celltype"), all.x = T)

merged_pos.df$p0_25 <- ifelse(merged_pos.df$expm1.x < merged_pos.df$expm1.y[,"p25.25%"], TRUE, FALSE)

merged_pos.df$p25_50 <- ifelse(merged_pos.df$expm1.x > merged_pos.df$expm1.y[,"p25.25%"] & 
                                 merged_pos.df$expm1.x < merged_pos.df$expm1.y[,"p50.50%"], TRUE, FALSE)

merged_pos.df$p50_75 <- ifelse(merged_pos.df$expm1.x > merged_pos.df$expm1.y[,"p50.50%"] & 
                                 merged_pos.df$expm1.x < merged_pos.df$expm1.y[,"p75.75%"], TRUE, FALSE)

merged_pos.df$p75_100 <- ifelse(merged_pos.df$expm1.x > merged_pos.df$expm1.y[,"p75.75%"], TRUE, FALSE)

merged_pos.df <- merged_pos.df[c("Cell_ID", "has.pSyn", "p0_25", "p25_50", "p50_75", "p75_100")]


# percentiles for pSyn- cells

pSyn_neg.df <- df[df$has.pSyn %in% "FALSE", ]

percentiles_neg <- aggregate(expm1 ~ orig.ident + celltype, data = pSyn_neg.df, FUN = function(x) {
  c(p25 = quantile(x, probs = 0.25),
    p50 = quantile(x, probs = 0.50),
    p75 = quantile(x, probs = 0.75))
})

merged_neg.df <- merge(pSyn_neg.df, percentiles_neg, by = c("orig.ident", "celltype"))

merged_neg.df$p0_25 <- ifelse(merged_neg.df$expm1.x < merged_neg.df$expm1.y[,"p25.25%"], TRUE, FALSE)

merged_neg.df$p25_50 <- ifelse(merged_neg.df$expm1.x > merged_neg.df$expm1.y[,"p25.25%"] & 
                                 merged_neg.df$expm1.x < merged_neg.df$expm1.y[,"p50.50%"], TRUE, FALSE)

merged_neg.df$p50_75 <- ifelse(merged_neg.df$expm1.x > merged_neg.df$expm1.y[,"p50.50%"] & 
                                 merged_neg.df$expm1.x < merged_neg.df$expm1.y[,"p75.75%"], TRUE, FALSE)

merged_neg.df$p75_100 <- ifelse(merged_neg.df$expm1.x > merged_neg.df$expm1.y[,"p75.75%"], TRUE, FALSE)

merged_neg.df <- merged_neg.df[c("Cell_ID", "has.pSyn", "p0_25", "p25_50", "p50_75", "p75_100")]

##

find_true_column <- function(row) {
  names(row)[which(row == TRUE)]
}

merged_pos.df$percentile <- apply(merged_pos.df[, c("p0_25", "p25_50", "p50_75", "p75_100")], 1, find_true_column)
merged_neg.df$percentile <- apply(merged_neg.df[, c("p0_25", "p25_50", "p50_75", "p75_100")], 1, find_true_column)

merged_all.df <- rbind(merged_pos.df, merged_neg.df)

for_xeniumexplorer.df <- merged_all.df[c("Cell_ID", "percentile")]
colnames(for_xeniumexplorer.df) <- c("cell_id", "group")
for_xeniumexplorer.df <- apply(for_xeniumexplorer.df, 2, as.character)

write.csv(for_xeniumexplorer.df, file = "/xenium_Line61_final/ctx/outs/tg_cells_percentiles_for_XE.csv")
