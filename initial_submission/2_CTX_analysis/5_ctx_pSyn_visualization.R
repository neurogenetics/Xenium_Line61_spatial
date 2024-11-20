library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
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

# for making cell labels within Xenium explorer (for Fig. 2a)
tg_1_metadata <- as.data.frame(xenium_ctx@meta.data)
write.csv(tg_1_metadata, file = "/xenium_Line61_final/ctx/outs/tg_1_metadata.csv")


# counting pSyn+/- by celltype

# table for pSyn counts by celltype by sample

sample_names <- xenium_ctx@meta.data$orig.ident
cell_types <- xenium_ctx@meta.data$celltype
pSyn_status <- xenium_ctx@meta.data$has_pSyn

pSyn_counts.df <- as.data.frame(table(sample_names, cell_types, pSyn_status))

pSyn_counts.df <- pSyn_counts.df[!grepl("^nTg", pSyn_counts.df$sample_names),, drop = FALSE]

total_counts <- pSyn_counts.df %>%
  group_by(sample_names, cell_types) %>%
  summarize(total_count = sum(Freq))

pSyn_counts.df <- merge(pSyn_counts.df, total_counts, by = c("sample_names", "cell_types"), all.x = TRUE)

pSyn_counts.df$pct <- pSyn_counts.df$Freq / pSyn_counts.df$total_count * 100

pSyn_counts.df <- pSyn_counts.df[!grepl("^FALSE", pSyn_counts.df$pSyn_status),, drop = FALSE]

selected_celltypes <- c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", "L5 ET ExN", "L5 NP ExN", "L6 IT ExN", "L6 CT ExN", "L6b ExN", 
                        "Pvalb+ InN", "Sst+ InN", "Lamp5+ InN", "Vip+ InN", 
                        "Micro", "Astro", "Oligo")

pSyn_counts.df <- pSyn_counts.df[pSyn_counts.df$cell_types %in% selected_celltypes, ]

pSyn_counts.df$cell_types <- factor(pSyn_counts.df$cell_types, levels = selected_celltypes)

saveRDS(pSyn_counts.df, file = "/xenium_Line61_final/ctx/outs/pSyn_counts.rds")
write.csv(pSyn_counts.df, file = "/xenium_Line61_final/ctx/outs/pSyn_counts.csv")

####################################################################################

pSyn_counts.df <- readRDS("/xenium_Line61_final/ctx/outs/pSyn_counts.rds")

####################################################################################

# using pSyn status as ident fov for rep image (for Fig. 2a)

# coords from selection in Xenium Explorer
crop <- Crop(xenium_ctx[["fov.5"]], x = c(4845.96, 5610.74), y = c(2397.31, 2893.61), coords = "tissue")
xenium_ctx[["crop"]] <- crop
DefaultBoundary(xenium_ctx[["crop"]]) <- "segmentation"

Idents(xenium_ctx) <- "has_pSyn"
levels(xenium_ctx) <- c("TRUE", "FALSE", "NA")
pSyn_colors <- c("red", "grey")

p2 <- ImageDimPlot(xenium_ctx, fov = "crop", dark.background = F, cols = pSyn_colors, axes = F, 
                   border.color = "black", border.size = 0.5) + NoLegend()

svglite("/xenium_Line61_final/ctx/outs/ctx_tg_1_zoom_pSyn_overlay.svg", 
        width = 5.9292, height = 8.8116)
plot(p2)
dev.off()

####################################################################################

# percentatges by cell type (Fig. 2b)

p3 <- ggplot(pSyn_counts.df, aes(x = cell_types, y = pct, fill = cell_types)) +
  geom_boxplot(outliers = F) +
  geom_jitter(width = 0.2) +
  scale_fill_manual(values = cluster_cols) +
  labs(x = " ",
       y = "% cells pSyn+") + 
  theme_bw() +
  guides(fill = "none") +
  theme(axis.text.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        text = element_text(family = "Arial")) +
  scale_x_discrete(limits = rev) + 
  coord_flip()

svglite("/xenium_Line61_final/ctx/outs/pct_pos_boxplot_flip.svg", 
        width = 7, height =8.25)
plot(p3)
dev.off()
