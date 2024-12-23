library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(extrafont)
library(svglite)

xenium_hip <- readRDS("/xenium_Line61_final/hip/xenium_HIP_metadata_pSyn.rds")


cluster_cols <- c("CA1" = "#98e86d", "CA2/3" = "#2fedcd", "DG" = "#ef7cf7", "InN" = "#9999f7", "Calb2+ ExN" = "#fae22d",
                  "Micro" = "#1a5402", "Astro" = "#2653ad", "Oligo" = "#999144", "OPC" = "#cc2400", "Endo" = "#ffa500", 
                  "VLMC" = "#f5259b", "VSMC" = "#7e14c9", "Cajal-Retzius" = "#9e9e9b")

####################################################################################

# overall umap reduction (Fig. 3b)

p1 <- DimPlot(xenium_hip, cols = cluster_cols, label = F, raster = F) + coord_fixed() + theme(text = element_text(family = "Arial"))

svglite("/xenium_Line61_final/hip/outs/hip_umap.svg", 
        width = 10, height = 8)
plot(p1)
dev.off()

####################################################################################

# cells visualized in space for 1 sample (Fig. 3a)

DefaultBoundary(xenium_hip[["fov"]]) <- "segmentation"

p2 <- ImageDimPlot(xenium_hip, fov = "fov", dark.background = F, cols = cluster_cols, axes = F, 
                   border.color = "black", border.size = 0.15) + NoLegend() + theme(text = element_text(family = "Arial"))

svglite("/xenium_Line61_final/hip/outs/hip_imagedimplot.svg", 
        width = 7.5, height = 8)
plot(p2)
dev.off()

####################################################################################

# celltype proportions plot (Fig. 3c)

library(speckle)
library(limma)

props.df <- as.data.frame(propeller(clusters = xenium_hip$celltype, sample = xenium_hip$orig.ident, group = xenium_hip$genotype, transform = "asin"))
props.df <- props.df[colnames(props.df) %in% c("BaselineProp.Freq", "BaselineProp.clusters")]
props.df <- props.df %>% 
  mutate(null = "1")

props.df$BaselineProp.clusters <- factor(props.df$BaselineProp.clusters, 
                                         levels = c("CA1", "CA2/3", "DG", "InN", "Calb2+ ExN", 
                                                    "Micro", "Astro", "Oligo", "OPC", "Endo", "VLMC",
                                                    "VSMC", "Cajal-Retzius"))


p3 <- ggplot(props.df, aes(x = null, y = BaselineProp.Freq, fill = BaselineProp.clusters)) + 
  geom_bar(position = position_stack(), stat = "identity", width = 0.5) + 
  labs(x = "",
       y = "Proportion of all cells") +
  scale_fill_manual(values = cluster_cols) +
  theme_bw() +
  ggtitle("HIP cell type proportions") + 
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

svglite("/xenium_Line61_final/hip/outs/hip_props.svg", 
        width = 4, height = 7)
plot(p3)
dev.off()

####################################################################################

# dotplot for marker genes (Fig. 3d)

levels(xenium_hip) <- rev(c("CA1", "CA2/3", "DG", "InN", "Calb2+ ExN", 
                            "Micro", "Astro", "Oligo", "OPC", "Endo", "VLMC",
                            "VSMC", "Cajal-Retzius"))

p4 <- DotPlot(xenium_hip, features = c("Slc17a7", "Pou3f1","Fibcd1","Slit2", "Prox1","Plekha2",
                                       "Gad2","Calb2","Siglech","Slc39a12","Opalin","Gpr17",
                                       "Cldn5", "Igf2", "Carmn", "Ebf3")) +
  RotatedAxis() + 
  theme(legend.position = "bottom", 
        legend.title.position = "top",
        text = element_text(family = "Arial")) +
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled), shape = 21, colour = "black", stroke = 0.5) +
  scale_fill_gradient2(low = "#05409e", mid = "white", high = "#e64a02") +
  guides(colour = "none",
         fill = guide_colourbar(title = "Avg. expression"),
         size = guide_legend(title = "Pct. expressed"))


svglite("/xenium_Line61_final/hip/outs/hip_dotplot.svg", 
        width = 8.5, height = 4.5)
plot(p4)
dev.off()

####################################################################################

hip_ctsbysample <- table(Idents(xenium_hip), xenium_hip$orig.ident)
write.csv(hip_ctsbysample,
          file = "/xenium_Line61_final/hip/outs/hip_ctsbysample.csv")
