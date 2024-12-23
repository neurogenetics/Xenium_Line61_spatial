library(ggplot2)
library(tidyverse)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(extrafont)


# load in general stuff

gene_order <- read.csv2("/xenium_Line61_final/xenium_upload1/general/gene_order_for_hmp.csv", sep = ",", header = F)

# general functions/settings

# create row split parameters for heatmap

row_split = rep("A_α-synuclein", 43)
row_split[3:6] = "B_α-syn phosphorylation"
row_split[7:9] = "C_Molecular chaperones"
row_split[10:17] = "D_Autophagy"
row_split[18:24] = "E_Ubiquitin-proteasome"
row_split[25:34] = "F_Endolysosome"
row_split[35:38] = "G_Mitochondria"
row_split[39:43] = "H_Other"


# set global heatmap options
ht_opt(RESET = TRUE)
ht_opt(legend_border = "black",
       heatmap_border = TRUE)


# set colors for heatmap

col_fun = colorRamp2(c(-0.5, 0, 0.5), c("#05409e", "white", "#e64a02"))
col_fun2 = colorRamp2(c(-8, 0, 4), c("#500c85", "white", "#ecf024"))

# set annotation for cell types

neurons <- c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", 
             "L5 ET ExN", "L5 NP ExN", "L6 IT ExN", 
             "L6 CT ExN", "L6b ExN")

cols <- c("L2/3 IT ExN" = "#fae22d", "L4/5 IT ExN" = "#44e3bb", 
          "L5 IT ExN" = "#80033d", "L5 ET ExN" = "#1e8a3b", 
          "L5 NP ExN" = "#82daf5", "L6 IT ExN" = "#c0ff4d", 
          "L6 CT ExN" = "#fa96fa", "L6b ExN" = "#278781")

ha = HeatmapAnnotation(`Cell type` = neurons,
                       col = list(`Cell type` = cols),
                       gp = gpar(col = "black"),
                       show_annotation_name = F,
                       annotation_legend_param = list(`Cell type` = list(nrow = 1,
                                                                         at = neurons)))

################################################################################

# 1. nTg vs. pSyn+ Tg

# make matrix for foldchanges

DE_for_hmp_pSyn <- readRDS("/xenium_Line61_final/ctx/outs/l2fc_nTg_Tg_pSyn_pos_VIZ.rds")
DE_for_hmp_pSyn.df <- as.data.frame(DE_for_hmp_pSyn)
DE_for_hmp_pSyn.df <- DE_for_hmp_pSyn.df[match(gene_order$V1, DE_for_hmp_pSyn.df$gene), ]

rownames(DE_for_hmp_pSyn.df) <- DE_for_hmp_pSyn.df[,1]
DE_for_hmp_pSyn.df[,1] <- NULL
rownames <- rownames(DE_for_hmp_pSyn.df)
DE_for_hmp_pSyn.df <- as.data.frame(apply(DE_for_hmp_pSyn.df, 2, FUN = as.numeric))
rownames(DE_for_hmp_pSyn.df) <- rownames

mat <- data.matrix(DE_for_hmp_pSyn.df)


# make matrix for FDRs

fdr_for_hmp_pSyn <- readRDS("/xenium_Line61_final/ctx/outs/fdr_nTg_Tg_pSyn_pos_VIZ.rds")
fdr_for_hmp_pSyn.df <- as.data.frame(fdr_for_hmp_pSyn)
fdr_for_hmp_pSyn.df <- fdr_for_hmp_pSyn.df[match(gene_order$V1, fdr_for_hmp_pSyn.df$gene), ]
fdr_for_hmp_pSyn.df[is.na(fdr_for_hmp_pSyn.df)] <- 1

rownames(fdr_for_hmp_pSyn.df) <- fdr_for_hmp_pSyn.df[,1]
fdr_for_hmp_pSyn.df[,1] <- NULL
rownames <- rownames(fdr_for_hmp_pSyn.df)
fdr_for_hmp_pSyn.df <- as.data.frame(apply(fdr_for_hmp_pSyn.df, 2, FUN = as.numeric))
rownames(fdr_for_hmp_pSyn.df) <- rownames

mat2 <- data.matrix(fdr_for_hmp_pSyn.df)


# set function for adding significance to hmp

fdr_fun1 <- function(j, i, x, y, w, h, fill) {
  if(mat2[i, j] <= 0.0001) {
    grid.text("﹡﹡﹡﹡", x, y)
  } else if(mat2[i, j] <= 0.001) {
    grid.text("﹡﹡﹡", x, y)
  } else if(mat2[i, j] <= 0.01) {
    grid.text("﹡﹡", x, y)
  } else if(mat2[i, j] <= 0.05) {
    grid.text("﹡", x, y)
  }
}


# make hmp

ht1 <- Heatmap(mat, 
               name = "log2FC", 
               col = col_fun, 
               cluster_rows = F, 
               cluster_columns = F, 
               row_names_side = "left",
               row_names_gp = gpar(fontsize = 12),
               row_split = row_split,
               row_title = c("", "α-syn phosphorylation", "Molecular chaperones",
                             "Autophagy", "Ubiquitin-proteasome", "Endolysosome", "Mitochondria", "Other"),
               row_title_gp = gpar(fontsize = 0),
               row_gap = unit(1.5, "mm"),
               row_title_rot = 0,
               row_title_side = "right", 
               heatmap_legend_param = list(direction = "horizontal", legend_width = unit(3, "cm"),
                                           title_position = "topleft"),
               cell_fun = fdr_fun1,
               column_title = "non-tg vs. pSyn+ α-syn-tg",
               column_title_gp = gpar(fontface = "bold"),
               column_names_gp = gpar(fontsize = 0),
               bottom_annotation = ha)

################################################################################

# 2. nTg vs. pSyn- Tg

# make matrix for foldchanges

DE_for_hmp_no_pSyn <- readRDS("/xenium_Line61_final/ctx/outs/l2fc_nTg_Tg_pSyn_neg_VIZ.rds")
DE_for_hmp_no_pSyn.df <- as.data.frame(DE_for_hmp_no_pSyn)
DE_for_hmp_no_pSyn.df <- DE_for_hmp_no_pSyn.df[match(gene_order$V1, DE_for_hmp_no_pSyn.df$gene), ]

rownames(DE_for_hmp_no_pSyn.df) <- DE_for_hmp_no_pSyn.df[,1]
DE_for_hmp_no_pSyn.df[,1] <- NULL
rownames <- rownames(DE_for_hmp_no_pSyn.df)
DE_for_hmp_no_pSyn.df <- as.data.frame(apply(DE_for_hmp_no_pSyn.df, 2, FUN = as.numeric))
rownames(DE_for_hmp_no_pSyn.df) <- rownames

mat3 <- data.matrix(DE_for_hmp_no_pSyn.df)


# make matrix for FDRs

fdr_for_hmp_no_pSyn <- readRDS("/xenium_Line61_final/ctx/outs/fdr_nTg_Tg_pSyn_neg_VIZ.rds")
fdr_for_hmp_no_pSyn.df <- as.data.frame(fdr_for_hmp_no_pSyn)
fdr_for_hmp_no_pSyn.df <- fdr_for_hmp_no_pSyn.df[match(gene_order$V1, fdr_for_hmp_no_pSyn.df$gene), ]
fdr_for_hmp_no_pSyn.df[is.na(fdr_for_hmp_no_pSyn.df)] <- 1

rownames(fdr_for_hmp_no_pSyn.df) <- fdr_for_hmp_no_pSyn.df[,1]
fdr_for_hmp_no_pSyn.df[,1] <- NULL
rownames <- rownames(fdr_for_hmp_no_pSyn.df)
fdr_for_hmp_no_pSyn.df <- as.data.frame(apply(fdr_for_hmp_no_pSyn.df, 2, FUN = as.numeric))
rownames(fdr_for_hmp_no_pSyn.df) <- rownames

mat4 <- data.matrix(fdr_for_hmp_no_pSyn.df)


# set function for adding significance to hmp

fdr_fun2 <- function(j, i, x, y, w, h, fill) {
  if(mat4[i, j] <= 0.0001) {
    grid.text("﹡﹡﹡﹡", x, y)
  } else if(mat4[i, j] <= 0.001) {
    grid.text("﹡﹡﹡", x, y)
  } else if(mat4[i, j] <= 0.01) {
    grid.text("﹡﹡", x, y)
  } else if(mat4[i, j] <= 0.05) {
    grid.text("﹡", x, y)
  }
}


# make hmp

ht2 <- Heatmap(mat3, 
               name = "log2FC", 
               col = col_fun, 
               cluster_rows = F, 
               cluster_columns = F, 
               row_names_side = "left",
               row_names_gp = gpar(fontsize = 12),
               row_split = row_split,
               row_title = c("", "α-syn phosphorylation", "Molecular chaperones",
                             "Autophagy", "Ubiquitin-proteasome", "Endolysosome", "Mitochondria", "Other"),
               row_title_gp = gpar(fontsize = 0),
               row_gap = unit(1.5, "mm"),
               row_title_rot = 0,
               row_title_side = "right", 
               heatmap_legend_param = list(direction = "horizontal", legend_width = unit(3, "cm"),
                                           title_position = "topleft"),
               cell_fun = fdr_fun2,
               column_title = "non-tg vs. pSyn- α-syn-tg",
               column_title_gp = gpar(fontface = "bold"),
               column_names_gp = gpar(fontsize = 0),
               bottom_annotation = ha)


################################################################################

# 3. Tg pSyn+ vs. pSyn-


# make matrix for foldchanges

DE_for_hmp_tg_only <- readRDS("/xenium_Line61_final/ctx/outs/l2fc_Tg_only_VIZ.rds")
DE_for_hmp_tg_only.df <- as.data.frame(DE_for_hmp_tg_only)
DE_for_hmp_tg_only.df <- DE_for_hmp_tg_only.df[match(gene_order$V1, DE_for_hmp_tg_only.df$gene), ]

rownames(DE_for_hmp_tg_only.df) <- DE_for_hmp_tg_only.df[,1]
DE_for_hmp_tg_only.df[,1] <- NULL
rownames <- rownames(DE_for_hmp_tg_only.df)
DE_for_hmp_tg_only.df <- as.data.frame(apply(DE_for_hmp_tg_only.df, 2, FUN = as.numeric))
rownames(DE_for_hmp_tg_only.df) <- rownames

mat5 <- data.matrix(DE_for_hmp_tg_only.df)


# make matrix for FDRs

fdr_for_hmp_tg_only <- readRDS("/xenium_Line61_final/ctx/outs/fdr_Tg_only_VIZ.rds")
fdr_for_hmp_tg_only.df <- as.data.frame(fdr_for_hmp_tg_only)
fdr_for_hmp_tg_only.df <- fdr_for_hmp_tg_only.df[match(gene_order$V1, fdr_for_hmp_tg_only.df$gene), ]
fdr_for_hmp_tg_only.df[is.na(fdr_for_hmp_tg_only.df)] <- 1

rownames(fdr_for_hmp_tg_only.df) <- fdr_for_hmp_tg_only.df[,1]
fdr_for_hmp_tg_only.df[,1] <- NULL
rownames <- rownames(fdr_for_hmp_tg_only.df)
fdr_for_hmp_tg_only.df <- as.data.frame(apply(fdr_for_hmp_tg_only.df, 2, FUN = as.numeric))
rownames(fdr_for_hmp_tg_only.df) <- rownames

mat6 <- data.matrix(fdr_for_hmp_tg_only.df)


# set function for adding significance to hmp

fdr_fun3 <- function(j, i, x, y, w, h, fill) {
  if(mat6[i, j] <= 0.0001) {
    grid.text("﹡﹡﹡﹡", x, y)
  } else if(mat6[i, j] <= 0.001) {
    grid.text("﹡﹡﹡", x, y)
  } else if(mat6[i, j] <= 0.01) {
    grid.text("﹡﹡", x, y)
  } else if(mat6[i, j] <= 0.05) {
    grid.text("﹡", x, y)
  }
}


# make hmp

ht3 <- Heatmap(mat5, 
               name = "log2FC", 
               col = col_fun, 
               cluster_rows = F, 
               cluster_columns = F, 
               row_names_side = "left",
               row_names_gp = gpar(fontsize = 12),
               row_split = row_split,
               row_title = c("", "α-syn phosphorylation", "Molecular chaperones",
                             "Autophagy", "Ubiquitin-proteasome", "Endolysosome", "Mitochondria", "Other"),
               row_title_gp = gpar(fontsize = 0),
               row_gap = unit(1.5, "mm"),
               row_title_rot = 0,
               row_title_side = "right", 
               heatmap_legend_param = list(direction = "horizontal", legend_width = unit(3, "cm"),
                                           title_position = "topleft"),
               cell_fun = fdr_fun3,
               column_title = "pSyn+ vs. pSyn- α-syn-tg",
               column_title_gp = gpar(fontface = "bold"),
               column_names_gp = gpar(fontsize = 0),
               bottom_annotation = ha)


################################################################################

# 4. GLM, hSNCA and gene expression

# make matrix for estimates

GLM_estimates <- readRDS("/xenium_Line61_final/ctx/outs/GLM_estimates.rds")
GLM_estimates.df <- as.data.frame(GLM_estimates)
GLM_estimates.df <- GLM_estimates.df[match(gene_order$V1, rownames(GLM_estimates.df)), ]

mat7 <- data.matrix(GLM_estimates.df)


# make matrix for FDRs

GLM_qvals <- readRDS("/xenium_Line61_final/ctx/outs/GLM_qvals.rds")
GLM_qvals.df <- as.data.frame(GLM_qvals)
GLM_qvals.df <- GLM_qvals.df[match(gene_order$V1, rownames(GLM_qvals.df)), ]

mat8 <- data.matrix(GLM_qvals.df)


# set function for adding significance to hmp

fdr_fun4 <- function(j, i, x, y, w, h, fill) {
  if(mat8[i, j] <= 0.0001) {
    grid.text("﹡﹡﹡﹡", x, y)
  } else if(mat8[i, j] <= 0.001) {
    grid.text("﹡﹡﹡", x, y)
  } else if(mat8[i, j] <= 0.01) {
    grid.text("﹡﹡", x, y)
  } else if(mat8[i, j] <= 0.05) {
    grid.text("﹡", x, y)
  }
}


# make hmp

ht4 <- Heatmap(mat7, 
               name = "GLM estimate (x10^4)", 
               col = col_fun2, 
               cluster_rows = F, 
               cluster_columns = F, 
               row_names_side = "left",
               row_names_gp = gpar(fontsize = 12),
               row_split = row_split,
               row_title = c("", "α-syn phosphorylation", "Molecular chaperones",
                             "Autophagy", "Ubiquitin-proteasome", "Endolysosome", "Mitochondria", "Other"),
               row_title_gp = gpar(fontsize = 0),
               row_gap = unit(1.5, "mm"),
               row_title_rot = 0,
               row_title_side = "right", 
               heatmap_legend_param = list(direction = "horizontal", legend_width = unit(3, "cm"),
                                           title_position = "topleft"),
               cell_fun = fdr_fun4,
               column_title = "GLM, GEx ~hSNCA",
               column_title_gp = gpar(fontface = "bold"),
               column_names_gp = gpar(fontsize = 0),
               bottom_annotation = ha)

##################################################################################################

# combine hmps 1+2+3+4 (Fig. 6b-e)

ht_list2 = ht1 + ht2 + ht3 + ht4

png("/xenium_Line61_final/ctx/outs/DE_hmp_allcomps_w_GLM.png", width = 12.8, height = 9.4, units = "in", res = 600)
draw(ht_list2, row_sub_title_side = "right", 
     annotation_legend_side = "bottom", heatmap_legend_side = "bottom", merge_legends = T)
dev.off()
