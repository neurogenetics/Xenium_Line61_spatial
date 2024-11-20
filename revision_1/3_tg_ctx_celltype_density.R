library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(extrafont)
library(svglite)
library(scCustomize)

xenium_ctx <- readRDS("/xenium_Line61_final/ctx/xenium_CTX_metadata_pSyn.rds")

cluster_cols <- c("L2/3 IT ExN" = "#fae22d", "L4/5 IT ExN" = "#44e3bb", "L5 IT ExN" = "#80033d", "L5 ET ExN" = "#1e8a3b", 
                  "L5 NP ExN" = "#82daf5", "L6 IT ExN" = "#c0ff4d", "L6 CT ExN" = "#fa96fa", "L6b ExN" = "#278781", 
                  "L2 IT RSP ExN" = "#3291fc","L4 RSP-ACA ExN" = "#c96004", "SUB-ProS ExN" = "#940099", "Car3 ExN" = "#925c99", 
                  "Pvalb+ InN" = "#99a2ff", "Sst+ InN" = "#ff9999", "Lamp5+ InN" = "#0f3466", "Vip+ InN" = "#967b09",
                  "Micro" = "#1a5402", "Astro" = "#2653ad", "Oligo" = "#999144", "OPC" = "#cc2400", "Endo" = "#ffa500", 
                  "VLMC" = "#f5365c", "VSMC" = "#b505eb", "Epen" = "#10e03d")

####################################################################################

crop <- Crop(xenium_ctx[["fov.5"]], x = c(4600, 5600), y = c(1800, 2600), coords = "tissue")
xenium_ctx[["crop"]] <- crop
DefaultBoundary(xenium_ctx[["crop"]]) <- "segmentation"

p1 <- ImageDimPlot(xenium_ctx, fov = "crop", dark.background = F, size = 1.5, cols = cluster_cols, axes = T, 
                   border.color = "black", border.size = 0.5) + 
  NoLegend() + 
  theme_bw() +
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 18)) +
  guides(fill = "none")

svglite("/xenium_Line61_final/ctx/outs/REVISED_TG_ctx_imagedimplot_zoom_w_axes.svg", 
        width = 9, height = 8)
plot(p1)
dev.off()

####################################################################################

# density plot for cell types by depth in the "crop" fov

centroids <- xenium_ctx@images$crop@boundaries$centroids
cells <- centroids@cells
coords <- as.data.frame(centroids@coords)

x_coords <- coords$x

df <- data.frame(cells = cells,
                 y = x_coords)

metadata <- xenium_ctx@meta.data
metadata <- metadata["celltype"]
metadata <- rownames_to_column(metadata, var = "cells")
metadata <- metadata[metadata$cells %in% df$cells, ]

df <- merge(df, metadata, by = "cells", all.x = T, all.y = T)

df <- df[df$celltype %in% c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", "L5 ET ExN", 
                            "L5 NP ExN", "L6 IT ExN", "L6 CT ExN", "L6b ExN", 
                            "Pvalb+ InN", "Sst+ InN", "Lamp5+ InN", "Vip+ InN", 
                            "Micro", "Astro", "Oligo", "OPC", "Endo", "VLMC", "VSMC", "Epen"), ]

df$celltype <- factor(df$celltype, levels = c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", "L5 ET ExN", 
                                              "L5 NP ExN", "L6 IT ExN", "L6 CT ExN", "L6b ExN", 
                                              "Pvalb+ InN", "Sst+ InN", "Lamp5+ InN", "Vip+ InN", 
                                              "Micro", "Astro", "Oligo", "OPC", "Endo", "VLMC", "VSMC", "Epen"))

p2 <- ggplot(df, aes(y = y)) + 
  geom_density(aes(colour = celltype, fill = celltype, alpha = 0.5)) +
  facet_wrap(~celltype, ncol = 4) +
  theme_bw() +
  scale_color_manual(values = cluster_cols) +
  scale_fill_manual(values = cluster_cols) + 
  labs(x = "Density",
       y = "Depth") +
  theme(text = element_text(family = "Arial"), 
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  guides(alpha = "none",
         colour = "none",
         fill = "none")

svglite("/xenium_Line61_final/ctx/outs/REVISION_TG_1ctx_densityplot.svg", 
        width = 6, height = 10)
plot(p2)
dev.off()
