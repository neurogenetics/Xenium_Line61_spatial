library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(broom)
library(ggrepel)
library(extrafont)
library(svglite)


xenium_hip <- readRDS("/xenium_Line61_final/hip/xenium_HIP_metadata_pSyn.rds")

cluster_cols <- c("CA1" = "#98e86d", "CA2/3" = "#2fedcd", "DG" = "#ef7cf7", "InN" = "#9999f7", "Calb2+ ExN" = "#fae22d",
                  "Micro" = "#1a5402", "Astro" = "#2653ad", "Oligo" = "#999144", "OPC" = "#cc2400", "Endo" = "#ffa500", 
                  "VLMC" = "#f5259b", "VSMC" = "#7e14c9", "Cajal-Retzius" = "#9e9e9b")

################################################################################

# generate average expression matrix w/ selected celltypes

Idents(xenium_hip) <- "celltype"
AvgExp_hip <- AverageExpression(xenium_hip, assays = "Xenium", group.by = c("ident", "orig.ident"), slot = "data")
AvgExp_hip.df <- as.data.frame(AvgExp_hip)

AvgExp_hip.df <- AvgExp_hip.df %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(cols = -gene, names_to = "Animal_CellType", values_to = "Expression")

AvgExp_hip.df <- AvgExp_hip.df %>%
  separate(Animal_CellType, into = c("celltype", "animal"), sep = "_", remove = T)

celltypes_keep <- c("Xenium.DG", "Xenium.CA2.3", "Xenium.CA1",
                    "Xenium.InN", "Xenium.Micro", "Xenium.Astro", 
                    "Xenium.Oligo")

AvgExp_hip.df <- AvgExp_hip.df[AvgExp_hip.df$celltype %in% celltypes_keep, ]

AvgExp_hip.df$celltype <- ifelse(AvgExp_hip.df$celltype == "Xenium.DG", "DG", 
                                 ifelse(AvgExp_hip.df$celltype == "Xenium.CA2.3", "CA2/3",
                                        ifelse(AvgExp_hip.df$celltype == "Xenium.CA1", "CA1",
                                               ifelse(AvgExp_hip.df$celltype == "Xenium.InN", "InN",
                                                      ifelse(AvgExp_hip.df$celltype == "Xenium.Micro", "Micro",
                                                             ifelse(AvgExp_hip.df$celltype == "Xenium.Astro", "Astro",
                                                                    ifelse(AvgExp_hip.df$celltype == "Xenium.Oligo", "Oligo", "Other")))))))

AvgExp_hip.df$celltype <- factor(AvgExp_hip.df$celltype, levels = c("CA1", "CA2/3", "DG", "InN", "Micro", "Astro", "Oligo")) 

saveRDS(AvgExp_hip.df, file = "/xenium_Line61_final/hip/outs/avgexp_hip_df.rds")

AvgExp_hip.df <- readRDS("/xenium_Line61_final/hip/outs/avgexp_hip_df.rds")

################################################################################

# hSNCA (Fig. 4c)

hsnca.df <- AvgExp_hip.df[AvgExp_hip.df$gene %in% "hSNCA", ]
hsnca.df$genotype <- ifelse(hsnca.df$animal == c("nTg.1", "nTg.2", "nTg.3", "nTg.4"), "non-tg", "α-syn-tg")

hsnca.df$celltype <- factor(hsnca.df$celltype, levels = c("CA1", "CA2/3", "DG", "InN", "Micro", "Astro", "Oligo")) 

p1 <- ggplot(hsnca.df, aes(x = genotype, y = Expression, fill = celltype)) +
  geom_boxplot(outliers = F, width = 0.9) + 
  geom_point(position = position_dodge(width = 0.9)) +
  labs(x = "Genotype", 
       y = "Average hSNCA expression \n (expm1[log-norm counts])") +
  ggtitle("hSNCA expression (by sample)") + 
  theme_bw() + 
  scale_fill_manual(values = cluster_cols,
                    labels = c("CA1", "CA2/3", "DG", "InN", "Micro", 
                               "Astro", "Oligo")) + 
  guides(fill = guide_legend(title = "Cell type")) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        text = element_text(family = "Arial")) +
  scale_y_continuous(limits = c(0, 700)) + 
  scale_x_discrete(expand = c(0,0))

svglite("/xenium_Line61_final/hip/outs/hip_hsnca_exp_by_sample_nTg_Tg.svg", 
        width = 9.5, height = 3.8)
plot(p1)
dev.off()

ggsave(p1, filename = "/xenium_Line61_final/hip/outs/hip_hsnca_exp_by_sample_nTg_Tg.png",
       width = 9.5, height = 3.8)


################################################################################

# Plk2 (Fig. 4d)

plk2.df <- AvgExp_hip.df[AvgExp_hip.df$gene %in% "Plk2", ]
plk2.df <- plk2.df[plk2.df$animal %in% c("nTg.1", "nTg.2", "nTg.3", "nTg.4"), ]

plk2.df$celltype <- factor(plk2.df$celltype, levels = c("CA1", "CA2/3", "DG", "InN", "Micro", "Astro", "Oligo")) 

p2 <- ggplot(plk2.df, aes(x = celltype, y = Expression, fill = celltype)) +
  geom_boxplot(outliers = F, width = 0.9) + 
  geom_point(position = position_dodge(width = 0.9)) +
  labs(x = "", 
       y = "Average Plk2 expression \n (expm1[log-norm counts])") +
  ggtitle("Plk2 expression (by sample, non-tg only)") + 
  theme_bw() + 
  scale_fill_manual(values = cluster_cols) + 
  guides(fill = "none") + 
  scale_x_discrete(labels = c("CA1", "CA2/3", "DG", "InN", "Micro", 
                              "Astro", "Oligo")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        text = element_text(family = "Arial")) +
  scale_y_continuous(limits = c(0, 250))  

svglite("/xenium_Line61_final/hip/outs/hip_plk2_exp_by_sample_nTg_Tg.svg", 
        width = 7.7, height = 3.3)
plot(p2)
dev.off()
