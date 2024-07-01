library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(broom)
library(ggrepel)
library(extrafont)
library(svglite)


xenium_ctx <- readRDS("/xenium_Line61_final/ctx/xenium_CTX_metadata_pSyn.rds")

cluster_cols <- c("L2/3 IT ExN" = "#fae22d", "L4/5 IT ExN" = "#44e3bb", "L5 IT ExN" = "#80033d", "L5 ET ExN" = "#1e8a3b", 
                  "L5 NP ExN" = "#82daf5", "L6 IT ExN" = "#c0ff4d", "L6 CT ExN" = "#fa96fa", "L6b ExN" = "#278781", 
                  "L2 IT RSP ExN" = "#3291fc","L4 RSP-ACA ExN" = "#c96004", "SUB-ProS ExN" = "#940099", "Car3 ExN" = "#925c99", 
                  "Pvalb+ InN" = "#99a2ff", "Sst+ InN" = "#ff9999", "Lamp5+ InN" = "#0f3466", "Vip+ InN" = "#967b09",
                  "Micro" = "#1a5402", "Astro" = "#2653ad", "Oligo" = "#999144", "OPC" = "#cc2400", "Endo" = "#ffa500", 
                  "VLMC" = "#f5365c", "VSMC" = "#b505eb", "Epen" = "#10e03d")

################################################################################

# generate average expression matrix w/ selected celltypes

Idents(xenium_ctx) <- "celltype"
AvgExp_ctx <- AverageExpression(xenium_ctx, assays = "Xenium", group.by = c("ident", "orig.ident"), slot = "data")
AvgExp_ctx.df <- as.data.frame(AvgExp_ctx)

AvgExp_ctx.df <- AvgExp_ctx.df %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(cols = -gene, names_to = "Animal_CellType", values_to = "Expression")

AvgExp_ctx.df <- AvgExp_ctx.df %>%
  separate(Animal_CellType, into = c("celltype", "animal"), sep = "_", remove = T)

celltypes_keep <- c("Xenium.L2.3.IT.ExN", "Xenium.L4.5.IT.ExN", "Xenium.L5.IT.ExN", "Xenium.L5.ET.ExN", "Xenium.L5.NP.ExN",
                    "Xenium.L6.IT.ExN", "Xenium.L6.CT.ExN", "Xenium.L6b.ExN", 
                    "Xenium.Pvalb..InN", "Xenium.Sst..InN", "Xenium.Lamp5..InN", "Xenium.Vip..InN",
                    "Xenium.Micro", "Xenium.Astro", "Xenium.Oligo")

AvgExp_ctx.df <- AvgExp_ctx.df[AvgExp_ctx.df$celltype %in% celltypes_keep, ]

renamed_celltypes <- c("Xenium.L2.3.IT.ExN" = "L2/3 IT ExN",
                       "Xenium.L4.5.IT.ExN" = "L4/5 IT ExN",
                       "Xenium.L5.IT.ExN" = "L5 IT ExN",
                       "Xenium.L5.ET.ExN" = "L5 ET ExN",
                       "Xenium.L5.NP.ExN" = "L5 NP ExN",
                       "Xenium.L6.IT.ExN" = "L6 IT ExN",
                       "Xenium.L6.CT.ExN" = "L6 CT ExN",
                       "Xenium.L6b.ExN" = "L6b ExN",
                       "Xenium.Pvalb..InN" = "Pvalb+ InN",
                       "Xenium.Sst..InN" = "Sst+ InN",
                       "Xenium.Lamp5..InN" = "Lamp5+ InN",
                       "Xenium.Vip..InN" = "Vip+ InN",
                       "Xenium.Micro" = "Micro",
                       "Xenium.Astro" = "Astro",
                       "Xenium.Oligo" = "Oligo")

AvgExp_ctx.df <- AvgExp_ctx.df %>%
  mutate(celltype = recode(celltype, !!!renamed_celltypes))

AvgExp_ctx.df$celltype <- factor(AvgExp_ctx.df$celltype, levels = c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", "L5 ET ExN", "L5 NP ExN", 
                                                                    "L6 IT ExN", "L6 CT ExN", "L6b ExN", 
                                                                    "Pvalb+ InN", "Sst+ InN", "Lamp5+ InN", "Vip+ InN", 
                                                                    "Micro", "Astro", "Oligo")) 

saveRDS(AvgExp_ctx.df, file = "/xenium_Line61_final/ctx/outs/avgexp_ctx_df.rds")

################################################################################

AvgExp_ctx.df <- readRDS("/xenium_Line61_final/ctx/outs/avgexp_ctx_df.rds")

################################################################################

# hSNCA expression by celltype (nTg and Tg) (Fig. 2d)

hsnca.df <- AvgExp_ctx.df[AvgExp_ctx.df$gene %in% "hSNCA", ]
hsnca.df$genotype <- ifelse(hsnca.df$animal == c("nTg.1", "nTg.2", "nTg.3", "nTg.4"), "non-tg", "α-syn-tg")

hsnca.df$celltype <- factor(hsnca.df$celltype, levels = c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", "L5 ET ExN", "L5 NP ExN", 
                                                          "L6 IT ExN", "L6 CT ExN", "L6b ExN", 
                                                          "Pvalb+ InN", "Sst+ InN", "Lamp5+ InN", "Vip+ InN", 
                                                          "Micro", "Astro", "Oligo")) 

p1 <- ggplot(hsnca.df, aes(x = genotype, y = Expression, fill = celltype)) +
  geom_boxplot(outliers = F, width = 0.9) + 
  geom_point(position = position_dodge(width = 0.9)) +
  labs(x = "Genotype", 
       y = "Average hSNCA expression \n (expm1[log-norm counts])") +
  ggtitle("hSNCA expression (by sample)") + 
  theme_bw() + 
  scale_fill_manual(values = cluster_cols,
                    labels = c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", "L5 ET ExN", "L5 NP ExN", 
                               "L6 IT ExN", "L6 CT ExN", "L6b ExN", 
                               "Pvalb+ InN", "Sst+ InN", "Lamp5+ InN", "Vip+ InN", 
                               "Micro", "Astro", "Oligo")) + 
  guides(fill = guide_legend(title = "Cell type")) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        text = element_text(family = "Arial")) +
  scale_y_continuous(limits = c(0, 750)) + 
  scale_x_discrete(expand = c(0,0))

svglite("/xenium_Line61_final/ctx/outs/ctx_hsnca_exp_by_sample_nTg_Tg.svg", 
        width = 14, height = 4)
plot(p1)
dev.off()

################################################################################

# Plk2 (Fig. 2e)

plk2.df <- AvgExp_ctx.df[AvgExp_ctx.df$gene %in% "Plk2", ]
plk2.df <- plk2.df[plk2.df$animal %in% c("nTg.1", "nTg.2", "nTg.3", "nTg.4"), ]

plk2.df$celltype <- factor(plk2.df$celltype, levels = c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", "L5 ET ExN", "L5 NP ExN", 
                                                        "L6 IT ExN", "L6 CT ExN", "L6b ExN", 
                                                        "Pvalb+ InN", "Sst+ InN", "Lamp5+ InN", "Vip+ InN", 
                                                        "Micro", "Astro", "Oligo")) 

p2 <- ggplot(plk2.df, aes(x = celltype, y = Expression, fill = celltype)) +
  geom_boxplot(outliers = F, width = 0.9) + 
  geom_point(position = position_dodge(width = 0.9)) +
  labs(x = "", 
       y = "Average Plk2 expression \n (expm1[log-norm counts])") +
  ggtitle("Plk2 expression (by sample, non-tg only)") + 
  theme_bw() + 
  scale_fill_manual(values = cluster_cols) + 
  guides(fill = "none") + 
  scale_x_discrete(labels = c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", "L5 ET ExN", "L5 NP ExN", 
                              "L6 IT ExN", "L6 CT ExN", "L6b ExN", 
                              "Pvalb+ InN", "Sst+ InN", "Lamp5+ InN", "Vip+ InN", 
                              "Micro", "Astro", "Oligo")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold"),
        text = element_text(family = "Arial")) +
  scale_y_continuous(limits = c(0, 150))  

svglite("/xenium_Line61_final/ctx/outs/ctx_plk2_exp_by_sample_nTg_Tg.svg", 
        width = 14, height = 4.5)
plot(p2)
dev.off()

################################################################################

# hSNCA featureplots b/w non-tg and tg (Fig. 2c)

p3 <- ImageFeaturePlot(xenium_ctx, features = "hSNCA", fov = c("fov", "fov.5"), 
                       dark.background = F, scale = "all", coord.fixed = T, combine = F) 
p4 <- p3[[1]] + ggtitle("non-tg") + NoLegend()
p5 <- p3[[2]] + ggtitle("α-syn-tg")
p6 <- p4 + p5 + theme(text = element_text(family = "Arial"))

svglite("/xenium_Line61_final/ctx/outs/hsnca_imagefeatplot.svg", 
        width = 9, height = 3)
plot(p6)
dev.off()
