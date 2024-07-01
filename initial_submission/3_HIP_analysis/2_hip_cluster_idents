library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)

xenium_hip <- readRDS("/xenium_Line61_final/hip/xenium_HIP-ONLY_lognorm_integrated_clustered.rds")

####################################################################################

# running FindAllMarkers

markers <- FindAllMarkers(xenium_hip, logfc.threshold = 0.5, only.pos = T, min.pct = 0.7)
write.csv(markers, file = "/xenium_Line61_final/hip/outs/xenium_hip_markers.csv")

####################################################################################

# renaming/releveling idents

Idents(xenium_hip) <- "hip.clusters"

xenium_hip <- RenameIdents(xenium_hip, `0` = "DG", `1` = "CA2/3", `2` = "CA1", `3` = "Astro", `4` = "Endo", `5` = "InN", 
                           `6` = "Oligo", `7` = "Micro", `8` = "OPC", `9` = "VLMC", `10` = "VSMC", 
                           `11` = "Cajal-Retzius", `12` = "Calb2+ ExN")


xenium_hip$celltype <- Idents(xenium_hip)
Idents(xenium_hip) <- "celltype"

Idents(xenium_hip) <- factor(x = Idents(xenium_hip), 
                             levels = c("CA1", "CA2/3", "DG", "InN", "Calb2+ ExN", "Micro", "Astro", 
                                        "Oligo", "OPC", "Endo", "VLMC", "VSMC", "Cajal-Retzius"))

####################################################################################

# set colors

cluster_colors <- c("CA1" = "#98e86d", "CA2/3" = "#2fedcd", "DG" = "#ef7cf7", "InN" = "#9999f7", "Calb2+ ExN" = "#fae22d",
                    "Micro" = "#1a5402", "Astro" = "#2653ad", "Oligo" = "#999144", "OPC" = "#cc2400", "Endo" = "#ffa500", 
                    "VLMC" = "#f5259b", "VSMC" = "#7e14c9", "Cajal-Retzius" = "#9e9e9b")


cluster_cols <- sapply(xenium_hip@meta.data$celltype, function(celltype) {
  cluster_colors[[as.character(celltype)]]
})

xenium_hip$cluster_cols <- cluster_cols

####################################################################################

# add metadata for genotype (nTg or Tg)

genotype <- ifelse(grepl("_[1-4]$", colnames(xenium_hip)), "nTg", "Tg")
xenium_hip <- AddMetaData(object = xenium_hip, metadata = genotype, col.name = "genotype")


# add metadata for sample (orig.ident)

orig.ident <- ifelse(grepl("_1", colnames(xenium_hip)), "nTg_1",
                     ifelse(grepl("_2", colnames(xenium_hip)), "nTg_2",
                            ifelse(grepl("_3", colnames(xenium_hip)), "nTg_3",
                                   ifelse(grepl("_4", colnames(xenium_hip)), "nTg_4",
                                          ifelse(grepl("_5", colnames(xenium_hip)), "Tg_1",
                                                 ifelse(grepl("_6", colnames(xenium_hip)), "Tg_2",
                                                        ifelse(grepl("_7", colnames(xenium_hip)), "Tg_3", "Tg_4")))))))
xenium_hip <- AddMetaData(object = xenium_hip, metadata = orig.ident, col.name = "orig.ident")

####################################################################################

saveRDS(xenium_hip, file = "/xenium_Line61_final/hip/xenium_HIP_metadata.rds")
