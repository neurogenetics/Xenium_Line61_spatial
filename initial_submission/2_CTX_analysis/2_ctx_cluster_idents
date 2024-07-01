library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)

xenium_ctx <- readRDS("/xenium_Line61_final/ctx/xenium_CTX-ONLY_lognorm_integrated_clustered.rds")

####################################################################################

# running FindAllMarkers

markers <- FindAllMarkers(xenium_ctx, logfc.threshold = 1, only.pos = T)
write.csv(markers, file = "/xenium_Line61_final/ctx/outs/ctx_markers.csv")


# imagedimplot highlighting 1 cell type for spatial identification

cols <- c("grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey",
          "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", 
          "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", 
          "grey", "grey", "grey", "red")

ImageDimPlot(xenium_ctx, fov = "fov", cols = cols, dark.background = F)

DimPlot(xenium_ctx, label = T)


# based on markers, clusters 31 - 34 need to be filtered (don't match any canonical cell types, may be poor segmentation)

Idents(xenium_ctx) <- "ctx.clusters"
xenium_ctx <- subset(xenium_ctx, idents = c("0","1","2","3","4","5","6","7","8","9","10",
                                            "11","12","13","14","15","16","17","18","19","20",
                                            "21","22","23","24","25","26","27","28","29","30"))

####################################################################################

# renaming/releveling idents

Idents(xenium_ctx) <- "ctx.clusters"

xenium_ctx <- RenameIdents(xenium_ctx, `0` = "L2/3 IT ExN", `1` = "L4/5 IT ExN", `2` = "L6 CT ExN", `3` = "Endo", `4` = "Oligo", `5` = "Astro", 
                           `6` = "Oligo", `7` = "L2 IT RSP ExN", `8` = "L5 IT ExN", `9` = "L6 IT ExN", `10` = "Micro", 
                           `11` = "L5 ET ExN", `12` = "Pvalb+ InN", `13` = "OPC", `14` = "Astro", `15` = "VLMC", 
                           `16` = "L5 NP ExN", `17` = "Sst+ InN", `18` = "VLMC", `19` = "L2/3 IT ExN", `20` = "VSMC", 
                           `21` = "Lamp5+ InN", `22` = "Vip+ InN", `23` = "Oligo", `24` = "Car3 ExN", `25` = "L6b ExN", 
                           `26` = "L4 RSP-ACA ExN", `27` = "SUB-ProS ExN", `28` = "L6 CT ExN", `29` = "Oligo", `30` = "Epen")


# 23 final celltypes after combining

xenium_ctx$celltype <- Idents(xenium_ctx)
Idents(xenium_ctx) <- "celltype"

Idents(xenium_ctx) <- factor(x = Idents(xenium_ctx), 
                             levels = c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", "L5 ET ExN", "L5 NP ExN", "L6 IT ExN", "L6 CT ExN", "L6b ExN", 
                                        "L2 IT RSP ExN", "L4 RSP-ACA ExN", "SUB-ProS ExN", "Car3 ExN",
                                        "Pvalb+ InN", "Sst+ InN", "Lamp5+ InN", "Vip+ InN",
                                        "Micro", "Astro", "Oligo", "OPC", "Endo", "VLMC", "VSMC", "Epen"))

####################################################################################

# set colors

cluster_colors <- c("L2/3 IT ExN" = "#fae22d", "L4/5 IT ExN" = "#44e3bb", "L5 IT ExN" = "#80033d", "L5 ET ExN" = "#1e8a3b", 
                    "L5 NP ExN" = "#82daf5", "L6 IT ExN" = "#c0ff4d", "L6 CT ExN" = "#fa96fa", "L6b ExN" = "#278781", 
                    "L2 IT RSP ExN" = "#3291fc","L4 RSP-ACA ExN" = "#c96004", "SUB-ProS ExN" = "#940099", "Car3 ExN" = "#925c99", 
                    "Pvalb+ InN" = "#99a2ff", "Sst+ InN" = "#ff9999", "Lamp5+ InN" = "#0f3466", "Vip+ InN" = "#967b09",
                    "Micro" = "#1a5402", "Astro" = "#2653ad", "Oligo" = "#999144", "OPC" = "#cc2400", "Endo" = "#ffa500", 
                    "VLMC" = "#f5365c", "VSMC" = "#b505eb", "Epen" = "#10e03d")


cluster_cols <- sapply(xenium_ctx@meta.data$celltype, function(celltype) {
  cluster_colors[[as.character(celltype)]]
})

xenium_ctx$cluster_cols <- cluster_cols

####################################################################################

# add metadata for genotype (nTg or Tg)

genotype <- ifelse(grepl("_[1-4]$", colnames(xenium_ctx)), "nTg", "Tg")
xenium_ctx <- AddMetaData(object = xenium_ctx, metadata = genotype, col.name = "genotype")


# add metadata for sample (orig.ident)

orig.ident <- ifelse(grepl("_1", colnames(xenium_ctx)), "nTg_1",
                     ifelse(grepl("_2", colnames(xenium_ctx)), "nTg_2",
                            ifelse(grepl("_3", colnames(xenium_ctx)), "nTg_3",
                                   ifelse(grepl("_4", colnames(xenium_ctx)), "nTg_4",
                                          ifelse(grepl("_5", colnames(xenium_ctx)), "Tg_1",
                                                 ifelse(grepl("_6", colnames(xenium_ctx)), "Tg_2",
                                                        ifelse(grepl("_7", colnames(xenium_ctx)), "Tg_3", "Tg_4")))))))
xenium_ctx <- AddMetaData(object = xenium_ctx, metadata = orig.ident, col.name = "orig.ident")

####################################################################################

saveRDS(xenium_ctx, 
        file = "/xenium_Line61_final/ctx/xenium_CTX_metadata.rds")
