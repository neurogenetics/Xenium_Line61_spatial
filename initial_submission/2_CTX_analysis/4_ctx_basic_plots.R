library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
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

####################################################################################

# overall umap reduction (Fig. 1c)

p1 <- DimPlot(xenium_ctx, cols = cluster_cols, label = F, raster = F) + coord_fixed() + theme(text = element_text(family = "Arial"))

svglite("/xenium_Line61_final/ctx/outs/ctx_umap.svg", 
        width = 10, height = 8)
plot(p1)
dev.off()

####################################################################################

# cells visualized in space for 1 sample (Fig. 1a)

DefaultBoundary(xenium_ctx[["fov"]]) <- "segmentation"

p2 <- ImageDimPlot(xenium_ctx, fov = "fov", dark.background = F, cols = cluster_cols, axes = F, 
                   border.color = "black", border.size = 0.15) + NoLegend() + theme(text = element_text(family = "Arial"))

svglite("/xenium_Line61_final/ctx/outs/ctx_imagedimplot.svg", 
        width = 10, height = 8)
plot(p2)
dev.off()

####################################################################################

# zoom in on cells in space for 1 sample (Fig. 1b)

crop <- Crop(xenium_ctx[["fov"]], x = c(5500, 7000), y = c(3000, 4000), coords = "tissue")
xenium_ctx[["crop"]] <- crop
DefaultBoundary(xenium_ctx[["crop"]]) <- "segmentation"

p3 <- ImageDimPlot(xenium_ctx, fov = "crop", dark.background = F, size = 1.5, cols = cluster_cols, axes = F, 
                   border.color = "black", border.size = 0.5) + NoLegend() + theme(text = element_text(family = "Arial"))

svglite("/xenium_Line61_final/ctx/outs/ctx_imagedimplot_zoom.svg", 
        width = 7, height = 6)
plot(p3)
dev.off()

####################################################################################

# celltype proportions plot (Fig. 1d)

library(speckle)
library(limma)

props.df <- as.data.frame(propeller(clusters = xenium_ctx$celltype, sample = xenium_ctx$orig.ident, group = xenium_ctx$genotype, transform = "asin"))
props.df <- props.df[colnames(props.df) %in% c("BaselineProp.Freq", "BaselineProp.clusters")]
props.df <- props.df %>% 
  mutate(null = "1")

props.df$BaselineProp.clusters <- factor(props.df$BaselineProp.clusters, 
                                         levels = c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", "L5 ET ExN", "L5 NP ExN", "L6 IT ExN", "L6 CT ExN", "L6b ExN", 
                                                    "L2 IT RSP ExN", "L4 RSP-ACA ExN", "SUB-ProS ExN", "Car3 ExN",
                                                    "Pvalb+ InN", "Sst+ InN", "Lamp5+ InN", "Vip+ InN",
                                                    "Micro", "Astro", "Oligo", "OPC", "Endo", "VLMC", "VSMC", "Epen"))


p4 <- ggplot(props.df, aes(x = null, y = BaselineProp.Freq, fill = BaselineProp.clusters)) + 
  geom_bar(position = position_stack(), stat = "identity", width = 0.5) + 
  labs(x = "",
       y = "Proportion of all cells") +
  scale_fill_manual(values = cluster_cols) +
  theme_bw() +
  ggtitle("CTX cell type proportions") + 
  theme(text = element_text(family ="Arial"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  guides(fill = guide_legend(ncol = 1)) + 
  scale_x_discrete(expand = expansion(mult = c(0, 0))) + 
  scale_y_continuous(expand = c(0,0))

svglite("/xenium_Line61_final/ctx/outs/ctx_props.svg", 
        width = 4, height = 7)
plot(p4)
dev.off()

####################################################################################

# dotplot for marker genes (Fig. 1e)

levels(xenium_ctx) <- rev(c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", "L5 ET ExN", "L5 NP ExN", "L6 IT ExN", "L6 CT ExN", "L6b ExN", 
                            "L2 IT RSP ExN", "L4 RSP-ACA ExN", "SUB-ProS ExN", "Car3 ExN",
                            "Pvalb+ InN", "Sst+ InN", "Lamp5+ InN", "Vip+ InN",
                            "Micro", "Astro", "Oligo", "OPC", "Endo", "VLMC", "VSMC", "Epen"))

p5 <- DotPlot(xenium_ctx, features = c("Slc17a7", "Lamp5", "Cux2", 
                                       "Kcnh5", "Rorb", 
                                       "Deptor", 
                                       "Bcl11b", "Gm19410", 
                                       "Vwc2l", "Nxph3",
                                       "Sema3e", 
                                       "Foxp2", "Rprm", 
                                       "Ccn2", "Cplx3", 
                                       "Rnf152",
                                       "Cbln1",
                                       "Nts",
                                       "Tmem163",
                                       "Gad1", "Pvalb", 
                                       "Sst", 
                                       "Vip", 
                                       "Laptm5", "Siglech",
                                       "Slc39a12", "Aqp4", 
                                       "Opalin", "Gjc3", 
                                       "Pdgfra", "Gpr17", 
                                       "Cldn5", "Ly6a", 
                                       "Igf2", "Fmod", 
                                       "Carmn", "Cspg4",
                                       "Spag16")) +
  RotatedAxis() + 
  theme(legend.position = "bottom", 
        legend.title.position = "top",
        text = element_text(family = "Arial")) +
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled), shape = 21, colour = "black", stroke = 0.5) +
  scale_fill_gradient2(low = "#05409e", mid = "white", high = "#e64a02") +
  guides(colour = "none",
         fill = guide_colourbar(title = "Avg. expression"),
         size = guide_legend(title = "Pct. expressed"))

svglite("/xenium_Line61_final/ctx/outs/ctx_dotplot.svg", 
        width = 12.5, height = 6)
plot(p5)
dev.off()

####################################################################################

# density plot for cell types by depth in the "crop" fov (Supplementary Fig. 4b)

crop <- Crop(xenium_ctx[["fov"]], x = c(5500, 7000), y = c(3000, 4000), coords = "tissue")
xenium_ctx[["crop"]] <- crop

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

cluster_cols <- c("L2/3 IT ExN" = "#fae22d", "L4/5 IT ExN" = "#44e3bb", "L5 IT ExN" = "#80033d", "L5 ET ExN" = "#1e8a3b", 
                  "L5 NP ExN" = "#82daf5", "L6 IT ExN" = "#c0ff4d", "L6 CT ExN" = "#fa96fa", "L6b ExN" = "#278781", 
                  "Pvalb+ InN" = "#99a2ff", "Sst+ InN" = "#ff9999", "Lamp5+ InN" = "#0f3466", "Vip+ InN" = "#967b09",
                  "Micro" = "#1a5402", "Astro" = "#2653ad", "Oligo" = "#999144", "OPC" = "#cc2400", "Endo" = "#ffa500", 
                  "VLMC" = "#f5365c", "VSMC" = "#b505eb", "Epen" = "#10e03d")


p6 <- ggplot(df, aes(y = y)) + 
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


svglite("/xenium_Line61_final/ctx/outs/ctx_densityplot.svg", 
        width = 6, height = 10)
plot(p6)
dev.off()

####################################################################################

# cropped image dimplot but w/ axes (Supplementary Fig. 4a)

crop <- Crop(xenium_ctx[["fov"]], x = c(5500, 7000), y = c(3000, 4000), coords = "tissue")
xenium_ctx[["crop"]] <- crop
DefaultBoundary(xenium_ctx[["crop"]]) <- "segmentation"

p7 <- ImageDimPlot(xenium_ctx, fov = "crop", dark.background = F, size = 1.5, cols = cluster_cols, axes = T, 
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

svglite("/xenium_Line61_final/ctx/outs/ctx_imagedimplot_zoom_w_axes.svg", 
        width = 9, height = 8)
plot(p7)
dev.off()

####################################################################################

ctx_ctsbysample <- table(Idents(xenium_ctx), xenium_ctx$orig.ident)
write.csv(ctx_ctsbysample,
          file = "/xenium_Line61_final/ctx/outs/ctx_ctsbysample.csv")
