library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(extrafont)
library(svglite)


xenium_hip <- readRDS("/xenium_Line61_final/hip/xenium_HIP_metadata_pSyn.rds")

cluster_cols <- c("CA1" = "#98e86d", "CA2/3" = "#2fedcd", "DG" = "#ef7cf7", "InN" = "#9999f7", "Calb2+ ExN" = "#fae22d",
                  "Micro" = "#1a5402", "Astro" = "#2653ad", "Oligo" = "#999144", "OPC" = "#cc2400", "Endo" = "#ffa500", 
                  "VLMC" = "#f5259b", "VSMC" = "#7e14c9", "Cajal-Retzius" = "#9e9e9b")

####################################################################################

# make 3 overlays in XE: CA1, CA2/3, DG (for Fig. 3a)
# for making cell labels within Xenium explorer
tg_1_metadata <- as.data.frame(xenium_hip@meta.data)
tg_1_metadata <- tg_1_metadata[tg_1_metadata$orig.ident %in% "Tg_1", ]
write.csv(tg_1_metadata, file = "/xenium_Line61_final/hip/outs/tg_1_metadata.csv")

####################################################################################

# counting pSyn+/- by celltype

# table for pSyn counts by celltype by sample

sample_names <- xenium_hip@meta.data$orig.ident
cell_types <- xenium_hip@meta.data$celltype
pSyn_status <- xenium_hip@meta.data$has_pSyn

pSyn_counts.df <- as.data.frame(table(sample_names, cell_types, pSyn_status))

pSyn_counts.df <- pSyn_counts.df[!grepl("^nTg", pSyn_counts.df$sample_names),, drop = FALSE]

total_counts <- pSyn_counts.df %>%
  group_by(sample_names, cell_types) %>%
  summarize(total_count = sum(Freq))

pSyn_counts.df <- merge(pSyn_counts.df, total_counts, by = c("sample_names", "cell_types"), all.x = TRUE)

pSyn_counts.df$pct <- pSyn_counts.df$Freq / pSyn_counts.df$total_count * 100

pSyn_counts.df <- pSyn_counts.df[!grepl("^FALSE", pSyn_counts.df$pSyn_status),, drop = FALSE]

selected_celltypes <- c("CA1", "CA2/3", "DG", "InN", "Micro", "Astro", "Oligo")

pSyn_counts.df <- pSyn_counts.df[pSyn_counts.df$cell_types %in% selected_celltypes, ]

pSyn_counts.df$cell_types <- factor(pSyn_counts.df$cell_types, levels = selected_celltypes)

saveRDS(pSyn_counts.df, file = "/xenium_Line61_final/hip/outs/pSyn_counts.rds")
write.csv(pSyn_counts.df, file = "/xenium_Line61_final/hip/outs/pSyn_counts.csv")

pSyn_counts.df <- readRDS("/xenium_Line61_final/hip/outs/pSyn_counts.rds")

####################################################################################

# ggplot for percentatges (Fig. 4b)

p1 <- ggplot(pSyn_counts.df, aes(x = cell_types, y = pct, fill = cell_types)) +
  geom_boxplot(outliers = F) +
  geom_jitter(width = 0.2) +
  scale_fill_manual(values = cluster_cols) +
  labs(x = " ",
       y = "% cells pSyn+") + 
  theme_bw() +
  guides(fill = "none") +
  theme(axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(family = "Arial"))

svglite("/xenium_Line61_final/hip/outs/hip_tg_cells_pSyn_pcts.svg", 
        width = 8, height = 3.5)
plot(p1)
dev.off()
