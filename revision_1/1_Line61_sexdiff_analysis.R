library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(extrafont)
library(svglite)

cluster_cols <- c("L2/3 IT ExN" = "#fae22d", "L4/5 IT ExN" = "#44e3bb", "L5 IT ExN" = "#80033d", "L5 ET ExN" = "#1e8a3b", 
                  "L5 NP ExN" = "#82daf5", "L6 IT ExN" = "#c0ff4d", "L6 CT ExN" = "#fa96fa", "L6b ExN" = "#278781", 
                  "L2 IT RSP ExN" = "#3291fc","L4 RSP-ACA ExN" = "#c96004", "SUB-ProS ExN" = "#940099", "Car3 ExN" = "#925c99", 
                  "Pvalb+ InN" = "#99a2ff", "Sst+ InN" = "#ff9999", "Lamp5+ InN" = "#0f3466", "Vip+ InN" = "#967b09",
                  "Micro" = "#1a5402", "Astro" = "#2653ad", "Oligo" = "#999144", "OPC" = "#cc2400", "Endo" = "#ffa500", 
                  "VLMC" = "#f5365c", "VSMC" = "#b505eb", "Epen" = "#10e03d")

################################################################################

AvgExp_ctx.df <- readRDS("/xenium_Line61_final/ctx/outs/avgexp_ctx_df.rds")

################################################################################

# hSNCA expression in Tg, labeled for M/F (CTX)

hsnca_tg.df <- AvgExp_ctx.df[AvgExp_ctx.df$gene %in% "hSNCA", ]
hsnca_tg.df <- hsnca_tg.df[hsnca_tg.df$animal == c("Tg.1", "Tg.2", "Tg.3", "Tg.4"), ]

hsnca_tg.df$sex <- ifelse(hsnca_tg.df$animal == c("Tg.1", "Tg.2"), "Female", "Male")

hsnca_tg.df$celltype <- factor(hsnca_tg.df$celltype, levels = c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", "L5 ET ExN", "L5 NP ExN", 
                                                                "L6 IT ExN", "L6 CT ExN", "L6b ExN", 
                                                                "Pvalb+ InN", "Sst+ InN", "Lamp5+ InN", "Vip+ InN", 
                                                                "Micro", "Astro", "Oligo")) 


p1 <- ggplot(hsnca_tg.df, aes(x = celltype, y = Expression, fill = celltype)) +
  geom_boxplot(outliers = F, width = 0.9) + 
  geom_point(position = position_dodge(width = 0.9), aes(shape = sex), size = 3) +
  labs(x = "Genotype", 
       y = "Average hSNCA expression") +
  ggtitle("hSNCA expression (by sample)") + 
  theme_bw() + 
  scale_fill_manual(values = cluster_cols,
                    labels = c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", "L5 ET ExN", "L5 NP ExN", 
                               "L6 IT ExN", "L6 CT ExN", "L6b ExN", 
                               "Pvalb+ InN", "Sst+ InN", "Lamp5+ InN", "Vip+ InN", 
                               "Micro", "Astro", "Oligo")) + 
  scale_shape_manual(values = c(16, 17)) +
  guides(fill = "none") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        text = element_text(family = "Arial")) #+
  scale_y_continuous(limits = c(0, 750)) 

svglite("/xenium_Line61_final/ctx/outs/hsnca_exp_by_sex.svg", 
        width = 16.5, height = 7.5)
plot(p1)
dev.off()

################################################################################

pSyn_counts.df <- readRDS("/xenium_Line61_final/ctx/outs/pSyn_counts.rds")

################################################################################

# pSyn percents, labeled by M/F (CTX)

pSyn_counts_sex.df <- pSyn_counts.df

pSyn_counts_sex.df$sex <- ifelse(pSyn_counts_sex.df$sample_names == c("Tg_1"), "Female", 
                                 ifelse(pSyn_counts_sex.df$sample_names == c("Tg_2"), "Female", "Male"))


p2 <- ggplot(pSyn_counts_sex.df, aes(x = cell_types, y = pct, fill = cell_types)) +
  geom_boxplot(outliers = F, width = 0.9) + 
  geom_point(position = position_dodge(width = 0.9), aes(shape = sex), size = 3) +
  labs(x = "Genotype", 
       y = "Average hSNCA expression") +
  ggtitle("hSNCA expression (by sample)") + 
  theme_bw() + 
  scale_fill_manual(values = cluster_cols,
                    labels = c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", "L5 ET ExN", "L5 NP ExN", 
                               "L6 IT ExN", "L6 CT ExN", "L6b ExN", 
                               "Pvalb+ InN", "Sst+ InN", "Lamp5+ InN", "Vip+ InN", 
                               "Micro", "Astro", "Oligo")) + 
  scale_shape_manual(values = c(16, 17)) +
  guides(fill = "none") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        text = element_text(family = "Arial")) 

svglite("/xenium_Line61_final/ctx/outs/pSyn_pcts_by_sex.svg", 
        width = 16.5, height = 7.5)
plot(p2)
dev.off()

################################################################################
################################################################################
################################################################################

cluster_cols <- c("CA1" = "#98e86d", "CA2/3" = "#2fedcd", "DG" = "#ef7cf7", "InN" = "#9999f7", "Calb2+ ExN" = "#fae22d",
                  "Micro" = "#1a5402", "Astro" = "#2653ad", "Oligo" = "#999144", "OPC" = "#cc2400", "Endo" = "#ffa500", 
                  "VLMC" = "#f5259b", "VSMC" = "#7e14c9", "Cajal-Retzius" = "#9e9e9b")

AvgExp_hip.df <- readRDS("/xenium_Line61_final/hip/outs/avgexp_hip_df.rds")

################################################################################

# hSNCA expression in Tg, labeled for M/F (HIP)

hsnca_tg.df <- AvgExp_hip.df[AvgExp_hip.df$gene %in% "hSNCA", ]
hsnca_tg.df <- hsnca_tg.df[hsnca_tg.df$animal == c("Tg.1", "Tg.2", "Tg.3", "Tg.4"), ]

hsnca_tg.df$sex <- ifelse(hsnca_tg.df$animal == c("Tg.1", "Tg.2"), "Female", "Male")

hsnca_tg.df$celltype <- factor(hsnca_tg.df$celltype, levels = c("CA1", "CA2/3", "DG", "InN", "Micro", "Astro", "Oligo")) 


p3 <- ggplot(hsnca_tg.df, aes(x = celltype, y = Expression, fill = celltype)) +
  geom_boxplot(outliers = F, width = 0.9) + 
  geom_point(position = position_dodge(width = 0.9), aes(shape = sex), size = 3) +
  labs(x = "Genotype", 
       y = "Average hSNCA expression") +
  ggtitle("hSNCA expression (by sample)") + 
  theme_bw() + 
  scale_fill_manual(values = cluster_cols,
                    labels = c("CA1", "CA2/3", "DG", "InN", "Micro", "Astro", "Oligo")) + 
  scale_shape_manual(values = c(16, 17)) +
  guides(fill = "none") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        text = element_text(family = "Arial")) +
scale_y_continuous(limits = c(0, 700)) 

svglite("/xenium_Line61_final/hip/outs/hip_hsnca_exp_by_sex.svg", 
        width = 9.5, height = 3.8)
plot(p3)
dev.off()

################################################################################

pSyn_counts_hip.df <- readRDS("/xenium_Line61_final/hip/outs/pSyn_counts.rds")

################################################################################

pSyn_counts_sex.df <- pSyn_counts_hip.df

pSyn_counts_sex.df$sex <- ifelse(pSyn_counts_sex.df$sample_names == c("Tg_1"), "Female", 
                                 ifelse(pSyn_counts_sex.df$sample_names == c("Tg_2"), "Female", "Male"))


p4 <- ggplot(pSyn_counts_sex.df, aes(x = cell_types, y = pct, fill = cell_types)) +
  geom_boxplot(outliers = F, width = 0.9) + 
  geom_point(position = position_dodge(width = 0.9), aes(shape = sex), size = 3) +
  labs(x = "Genotype", 
       y = "Average hSNCA expression") +
  ggtitle("hSNCA expression (by sample)") + 
  theme_bw() + 
  scale_fill_manual(values = cluster_cols,
                    labels = c("CA1", "CA2/3", "DG", "InN", "Micro", "Astro", "Oligo")) + 
  scale_shape_manual(values = c(16, 17)) +
  guides(fill = "none") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        text = element_text(family = "Arial")) 

svglite("/xenium_Line61_final/hip/outs/hip_pSyn_pcts_by_sex.svg", 
        width = 9.5, height = 3.8)
plot(p4)
dev.off()
