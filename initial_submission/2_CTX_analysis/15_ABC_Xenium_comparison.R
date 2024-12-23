library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(broom)
library(ggrepel)
library(extrafont)
library(svglite)
library(magrittr)

xenium_ctx <- readRDS("/xenium_Line61_final/ctx/xenium_CTX_metadata_pSyn.rds")

ABC_isoctx <- readRDS("/new_10Xv2/isocortex_all_lognorm_filtered_seurat.rds")

xenium_genes.df <- as.data.frame(read.csv2("/xenium_Line61_final/xenium_upload1/general/xenium_genes.csv", 
                                           sep = ",", header = T))
all_genes <- xenium_genes.df$ensembl_all
base_genes <- all_genes[1:248]
custom_genes <- all_genes[249:347]


################################################################################

# run AverageExpression

avgexp_ABC <- AverageExpression(ABC_isoctx, assays = "RNA", group.by = c("celltype") , slot = "data")
avgexp_ABC.df <- as.data.frame(avgexp_ABC)
avgexp_ABC.df <- avgexp_ABC.df[rownames(avgexp_ABC.df) %in% all_genes, ]
saveRDS(avgexp_ABC.df, file = "/xenium_Line61_final/ctx/ABC_comp/avgexp_ABC.rds")


Idents(xenium_ctx) <- "genotype"
ntg_keep <- WhichCells(xenium_ctx, idents = "nTg")
# using subset_opt function found here: https://github.com/alikhuseynov/add-on_R/blob/develop/R/subset_obj_seurat_v2.R 
# because Seurat is bugged and won't let you subset if there aren't some cells from every FOV in the subset object
subset_opt <- function(
    object = NULL, 
    subset = NULL, 
    cells = NULL, 
    idents = NULL,
    features = NULL,
    Update.slots = TRUE,
    Update.object = TRUE,
    ...)
{
  
  if (Update.slots) { 
    message("Updating object slots..")
    object %<>% UpdateSlots()
  }
  
  message("Cloing object..")
  obj_subset <- object
  
  # sanity check - use only cell ids (no indices)
  if (all(is.integer(cells))) { 
    cells <- Cells(obj_subset)[cells]
  }
  
  if (!missing(subset) || !is.null(idents)) {
    message("Extracting cells matched to `subset` and/or `idents`")
  }
  
  if (class(obj_subset) == "FOV") {
    message("object class is `FOV` ")
    cells <- Cells(obj_subset)
  } else if (!class(obj_subset) == "FOV" && !missing(subset)) {
    subset <- enquo(arg = subset)
    # cells to keep in the object
    cells <-
      WhichCells(object = obj_subset, 
                 cells = cells,
                 idents = idents,
                 expression = subset,
                 return.null = TRUE, ...)
  } else if (!class(obj_subset) == "FOV" && !is.null(idents)) {
    cells <-
      WhichCells(object = obj_subset, 
                 cells = cells,
                 idents = idents,
                 return.null = TRUE, ...)
  } else if (is.null(cells)) {
    cells <- Cells(obj_subset)
  }
  
  # added support for object class `FOV`
  if (class(obj_subset) == "FOV") {
    message("Matching cells for object class `FOV`..")
    cells_check <- any(obj_subset %>% Cells %in% cells)
  } else { 
    # check if cells are present in all FOV
    message("Matching cells in FOVs..")
    cells_check <-
      lapply(Images(obj_subset) %>% seq, 
             function(i) { 
               any(obj_subset[[Images(obj_subset)[i]]][["centroids"]] %>% Cells %in% cells) 
             }) %>% unlist
  }
  
  if (all(cells_check)) { 
    message("Cell subsets are found in all FOVs!", "\n",
            "Subsetting object..")
    obj_subset %<>% base::subset(cells = cells, 
                                 idents = idents,
                                 features = features,
                                 ...)
    # subset FOVs
    message("Subsetting FOVs..")
    fovs <- 
      lapply(Images(obj_subset) %>% seq, function(i) {
        base::subset(x = obj_subset[[Images(obj_subset)[i]]],
                     cells = cells, 
                     idents = idents, 
                     features = features, 
                     ...)
      })
    # replace subsetted FOVs
    for (i in fovs %>% seq) { obj_subset[[Images(object)[i]]] <- fovs[[i]] }
    
  } else { 
    # if cells are present only in one or several FOVs:
    # subset FOVs
    fovs <- 
      lapply(Images(obj_subset) %>% seq, function(i) {
        if (any(obj_subset[[Images(obj_subset)[i]]][["centroids"]] %>% Cells %in% cells)) {
          message("Cell subsets are found only in FOV: ", "\n", Images(obj_subset)[i])
          message("Subsetting Centroids..")
          base::subset(x = obj_subset[[Images(obj_subset)[i]]],
                       cells = cells, 
                       idents = idents, 
                       features = features, 
                       ...)
        }
      })
    # remove FOVs with no matching cells
    message("Removing FOVs where cells are NOT found: ", "\n", 
            paste0(Images(object)[which(!cells_check == TRUE)], "\n"))
    # replace subsetted FOVs
    for (i in fovs %>% seq) { obj_subset[[Images(object)[i]]] <- fovs[[i]] } 
    
    # subset final object
    message("..subset final object")
    obj_subset %<>% 
      base::subset(cells = cells,
                   idents = idents,
                   features = features, 
                   ...)
  }
  
  if (Update.object && !class(obj_subset) == "FOV") { 
    message("Updating object..")
    obj_subset %<>% UpdateSeuratObject() }
  
  message("Object is ready!")
  return(obj_subset)
  
}
xenium_ntg_only <- subset_opt(xenium_ctx, cells = ntg_keep) 
avgexp_xenium <- AverageExpression(xenium_ntg_only, assays = "Xenium", group.by = c("celltype") , slot = "data")
avgexp_xenium.df <- as.data.frame(avgexp_xenium)
saveRDS(avgexp_xenium.df, file = "/xenium_Line61_final/ctx/ABC_comp/avgexp_xenium.rds")


################################################################################

avgexp_ABC.df <- readRDS("/xenium_Line61_final/ctx/ABC_comp/avgexp_ABC.rds")
avgexp_xenium.df <- readRDS("/xenium_Line61_final/ctx/ABC_comp/avgexp_xenium.rds")

################################################################################

# reorganize dataframes


# ABC

genenames <- xenium_genes.df[c("gene_all", "ensembl_all")]

avgexp_ABC.df$ensembl <- rownames(avgexp_ABC.df)

avgexp_ABC.df <- avgexp_ABC.df %>%
  pivot_longer(cols = -ensembl, names_to = "celltype", values_to = "expression")

avgexp_ABC.df <- merge(avgexp_ABC.df, genenames, by.x = 'ensembl', by.y = "ensembl_all", all.x = TRUE)

avgexp_ABC.df <- avgexp_ABC.df[, c('gene_all', 'celltype', 'expression')]

avgexp_ABC.df <- subset(avgexp_ABC.df, celltype != "RNA.Other")

colnames(avgexp_ABC.df) <- c("gene", "celltype", "expression")


# make percentiles

grouped_data_ABC <- split(avgexp_ABC.df, avgexp_ABC.df$gene)

zscore_ABC <- lapply(grouped_data_ABC, function(gene_data) {
  gene_data$z_score <- scale(gene_data$expression)
  return(gene_data)
})

avgexp_ABC.df <- do.call(rbind, zscore_ABC)
avgexp_ABC.df$z_score <- as.numeric(avgexp_ABC.df$z_score)
avgexp_ABC.df$percentile <- pnorm(avgexp_ABC.df$z_score)
rownames(avgexp_ABC.df) <- NULL

saveRDS(avgexp_ABC.df, "/xenium_Line61_final/ctx/ABC_comp/avg_ABC_w_percentiles.rds")

################################################################################

# Xenium

avgexp_xenium.df <- avgexp_xenium.df %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(cols = -gene, names_to = "celltype", values_to = "expression")


# filter for only selected celltypes

celltypes_keep <- c("Xenium.L2.3.IT.ExN", "Xenium.L4.5.IT.ExN", "Xenium.L5.IT.ExN",
                    "Xenium.L5.ET.ExN", "Xenium.L5.NP.ExN", "Xenium.L6.IT.ExN", 
                    "Xenium.L6.CT.ExN", "Xenium.L6b.ExN", "Xenium.Pvalb..InN",
                    "Xenium.Sst..InN", "Xenium.Lamp5..InN", "Xenium.Vip..InN", "Xenium.Micro",
                    "Xenium.Astro", "Xenium.Oligo", "Xenium.OPC", "Xenium.VSMC", 
                    "Xenium.VLMC", "Xenium.Endo")

avgexp_xenium.df <- avgexp_xenium.df[avgexp_xenium.df$celltype %in% celltypes_keep, ]


grouped_data_xenium <- split(avgexp_xenium.df, avgexp_xenium.df$gene)

zscore_xenium <- lapply(grouped_data_xenium, function(gene_data) {
  gene_data$z_score <- scale(gene_data$expression)
  return(gene_data)
})

avgexp_xenium.df <- do.call(rbind, zscore_xenium)
avgexp_xenium.df$z_score <- as.numeric(avgexp_xenium.df$z_score)
avgexp_xenium.df$percentile <- pnorm(avgexp_xenium.df$z_score)

saveRDS(avgexp_xenium.df, "/xenium_Line61_final/ctx/ABC_comp/avg_Xenium_w_percentiles.rds")


################################################################################

# making comparison b/w Xenium and ABC expression profiles for major cell types

avg_ABC.df <- readRDS("/xenium_Line61_final/ctx/ABC_comp/avg_ABC_w_percentiles.rds")
avg_xenium.df <- readRDS("/xenium_Line61_final/ctx/ABC_comp/avg_Xenium_w_percentiles.rds")

avg_ABC.df$celltype <- gsub("^RNA\\.", "", avg_ABC.df$celltype)
avg_xenium.df$celltype <- gsub("^Xenium\\.", "", avg_xenium.df$celltype)


merged.df <- merge(avg_ABC.df, avg_xenium.df, by = c("celltype", "gene"), all.x = T, all.y = T)

merged.df <- merged.df[c("celltype", "gene", "percentile.x", "percentile.y")]


colnames(merged.df) <- c("celltype", "gene", "percentile.ABC", "percentile.Xenium")

merged.df <- na.omit(merged.df) # removes 3 genes (hSNCA, Skp1a/Skp1, Lyz2)

saveRDS(merged.df, "/xenium_Line61_final/ctx/ABC_comp/merged_percentiles.rds")


# rearrange to have columns as "celltype.platform"

rearranged_merged.df <- merged.df %>%
  pivot_longer(
    cols = starts_with("percentile"), 
    names_to = "platform", 
    values_to = "expression"
  ) %>%
  mutate(platform = gsub("percentile.", "", platform)) %>%
  unite("celltype_platform", celltype, platform, sep = ".") %>%
  pivot_wider(
    names_from = celltype_platform, 
    values_from = expression
  )

rearranged_merged.df <- column_to_rownames(rearranged_merged.df, var = "gene")

saveRDS(rearranged_merged.df, "/xenium_Line61_final/ctx/ABC_comp/rearranged_merged_df.rds")


################################################################################

rearranged_merged.df <- readRDS("/xenium_Line61_final/ctx/ABC_comp/rearranged_merged_df.rds")


# filter for genes only in base panel

xenium_genes.df <- as.data.frame(read.csv2("/xenium_Line61_final/xenium_upload1/general/xenium_genes.csv", 
                                           sep = ",", header = T))
all_genes <- xenium_genes.df$gene_all
base_genes <- all_genes[1:248]
rearranged_merged.df <- rearranged_merged.df[rownames(rearranged_merged.df) %in% base_genes, ]


# compute Pearson coefficients b/w Xenium and ABC celltype expression profiles 

abc_columns <- colnames(rearranged_merged.df)[grepl(".ABC", colnames(rearranged_merged.df))]
xenium_columns <- colnames(rearranged_merged.df)[grepl(".Xenium", colnames(rearranged_merged.df))]


correlation.df <- matrix(NA, nrow = length(xenium_columns), ncol = length(abc_columns))
rownames(correlation.df) <- xenium_columns
colnames(correlation.df) <- abc_columns

for (i in seq_along(xenium_columns)) {
  for (j in seq_along(abc_columns)) {
    cor_test <- cor.test(rearranged_merged.df[[xenium_columns[i]]], rearranged_merged.df[[abc_columns[j]]], method = "pearson")
    correlation.df[i, j] <- cor_test$estimate  
  }
}

correlation.df <- as.data.frame(correlation.df)

rownames(correlation.df) <- c("Astro", "Endo", "L2/3 IT ExN", "L4/5 IT ExN", "L5 ET ExN", "L5 IT ExN", "L5 NP ExN", "L6 CT ExN",
                              "L6 IT ExN", "L6b ExN", "Lamp5+ InN", "Micro", "Oligo", "OPC", "Pvalb+ InN", "Sst+ InN", "Vip+ InN",
                              "VLMC", "VSMC")

colnames(correlation.df) <- c("Astro", "Endo", "L2/3 IT ExN", "L4/5 IT ExN", "L5 ET ExN", "L5 IT ExN", "L5 NP ExN", "L6 CT ExN",
                              "L6 IT ExN", "L6b ExN", "Lamp5+ InN", "Micro", "Oligo", "OPC", "Pvalb+ InN", "Sst+ InN", "Vip+ InN",
                              "VLMC", "VSMC")

order <- c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", "L5 ET ExN", "L5 NP ExN", "L6 IT ExN", "L6 CT ExN", "L6b ExN",
           "Pvalb+ InN", "Sst+ InN", "Vip+ InN", "Lamp5+ InN", "Micro", "Astro", "Oligo", "OPC", "Endo", "VLMC", "VSMC")

correlation.df <- correlation.df[match(order, rownames(correlation.df)), ]
correlation.df <- correlation.df[, match(order, colnames(correlation.df))]

saveRDS(correlation.df, "/xenium_Line61_final/ctx/ABC_comp/correlation_pearson.rds")

################################################################################

# make a heatmap to show correlations (Supplementary Fig. 3c)

library(ComplexHeatmap)
library(circlize)

# set global heatmap options
ht_opt(RESET = TRUE)
ht_opt(legend_border = "black",
       heatmap_border = TRUE)


# set colors for heatmap

col_fun = colorRamp2(c(-1, 0, 1), c("#05409e", "white", "#e64a02"))


# function to label cells w/ matching rowname and colname (i.e., matched pairs of cell types)

label_fun <- function(j, i, x, y, width, height, fill) {
  if (rownames(mat)[i] == colnames(mat)[j]) {
    grid.text(sprintf("%.3f", mat[i, j]), x, y, gp = gpar(fontsize = 11, col = "black"))
  }
}


# make data matrix

correlation.df <- readRDS("/xenium_Line61_final/ctx/ABC_comp/correlation_pearson.rds")

mat <- data.matrix(correlation.df)

ht1 <- Heatmap(mat, 
               name = "Pearson R", 
               col = col_fun, 
               cluster_rows = F, 
               cluster_columns = F, 
               row_names_side = "left",
               row_names_gp = gpar(fontsize = 11),
               column_names_side = "bottom",
               column_names_gp = gpar(fontsize = 11),
               row_title = "Xenium (non-tg only)",
               row_title_rot = 90,
               row_title_side = "left", 
               row_title_gp = gpar(fontface = "bold"),
               column_title = "ABC atlas (subclass)",
               column_title_rot = 0,
               column_title_side = "bottom",
               column_title_gp = gpar(fontface = "bold"),
               cell_fun = label_fun)

png("/xenium_Line61_final/ctx/ABC_comp/outs/pearson_correlation.png", width = 11.4, height = 6.1, units = "in", res = 600)
draw(ht1)
dev.off()

pdf("/xenium_Line61_final/ctx/ABC_comp/outs/pearson_correlation.pdf", width = 11.4, height = 6.1)
draw(ht1)
dev.off()

################################################################################
################################################################################
################################################################################

# correlation of expression of individual genes

merged.df <- readRDS("/xenium_Line61_final/ctx/ABC_comp/merged_percentiles.rds")

xenium_genes.df <- as.data.frame(read.csv2("/xenium_Line61_final/xenium_upload1/general/xenium_genes.csv", 
                                           ",", header = T))
all_genes <- xenium_genes.df$gene_all
base_genes <- all_genes[1:248]
custom_genes <- all_genes[249:347]


################################################################################

# for genes in the base panel only (Supplementary Fig. 3a)

merged_base.df <- merged.df[merged.df$gene %in% base_genes, ]

lm_result <- lm(percentile.ABC ~ percentile.Xenium, data = merged_base.df)
tidy_result <- tidy(lm_result)

r_squared <- summary(lm_result)$r.squared
p_value <- tidy_result$`p.value`[2]

merged_base.df$celltype <- factor(merged_base.df$celltype, levels = c("L2.3.IT.ExN", "L4.5.IT.ExN", "L5.IT.ExN", "L5.ET.ExN", "L5.NP.ExN", 
                                                                      "L6.IT.ExN", "L6.CT.ExN", "L6b.ExN", 
                                                                      "Pvalb..InN", "Sst..InN", "Lamp5..InN", "Vip..InN", 
                                                                      "Micro", "Astro", "Oligo", "OPC", "Endo", "VLMC", "VSMC")) 

celltype_colors <- c("L2.3.IT.ExN" = "#FFD700",
                     "L4.5.IT.ExN" = "#228B22",
                     "L5.IT.ExN" = "#00CED1",
                     "L5.ET.ExN" = "#A52A2A",
                     "L5.NP.ExN" = "#D2B48C",
                     "L6.IT.ExN" = "#9370DB",
                     "L6.CT.ExN" = "#FFA500",
                     "L6b.ExN" = "#A4B6BA",
                     "Pvalb..InN" = "#CC789D",
                     "Sst..InN" = "#87CEEB",
                     "Vip..InN" = "#20B2AA",
                     "Lamp5..InN" = "#F52536",
                     "Micro" = "#5225BA", 
                     "Astro" = "#c72a9a", 
                     "Oligo" = "#591A2B", 
                     "OPC" = "#C4AE6A",
                     "Endo" = "#000080", 
                     "VLMC" = "#FFFF00", 
                     "VSMC" = "#1E90FF")


p1 <- ggplot(merged_base.df, aes(x = percentile.ABC, y = percentile.Xenium)) + 
  geom_point(aes(color = celltype)) + 
  geom_smooth(method='lm') + 
  labs(x = "Scaled expression - ABC atlas",
       y = "Scaled expression - Xenium (non-tg only)",
       title = "Expression correlation, Xenium base panel genes") + 
  scale_color_manual(values = celltype_colors,
                     labels = c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", "L5 ET ExN", 
                                "L5 NP ExN", "L6 IT ExN", "L6 CT ExN", "L6b ExN",
                                "Pvalb+ InN", "Sst+ InN", "Lamp5+ InN", "Vip+ InN",
                                "Micro", "Astro", "Oligo", "OPC", "Endo", "VLMC", "VSMC")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        text = element_text(family = "Arial")) +
  guides(color = guide_legend(title="Cell type"))

svglite("/xenium_Line61_final/ctx/ABC_comp/outs/ctx_ABC_comp_base_genes.svg", 
        width = 8, height = 5.5)
plot(p1)
dev.off()

################################################################################

# for genes in the custom panel only (Supplementary Fig. 3b)

merged_custom.df <- merged.df[merged.df$gene %in% custom_genes, ]

lm_result <- lm(percentile.ABC ~ percentile.Xenium, data = merged_custom.df)
tidy_result <- tidy(lm_result)

r_squared <- summary(lm_result)$r.squared
p_value <- tidy_result$`p.value`[2]

merged_custom.df$celltype <- factor(merged_custom.df$celltype, levels = c("L2.3.IT.ExN", "L4.5.IT.ExN", "L5.IT.ExN", "L5.ET.ExN", "L5.NP.ExN", 
                                                                          "L6.IT.ExN", "L6.CT.ExN", "L6b.ExN", 
                                                                          "Pvalb..InN", "Sst..InN", "Lamp5..InN", "Vip..InN", 
                                                                          "Micro", "Astro", "Oligo", "OPC", "Endo", "VLMC", "VSMC")) 


p2 <- ggplot(merged_custom.df, aes(x = percentile.ABC, y = percentile.Xenium)) + 
  geom_point(aes(color = celltype)) + 
  geom_smooth(method='lm') + 
  labs(x = "Scaled expression - ABC atlas",
       y = "Scaled expression - Xenium (non-tg only)",
       title = "Expression correlation, Xenium custom panel genes") + 
  scale_color_manual(values = celltype_colors,
                     labels = c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", "L5 ET ExN", 
                                "L5 NP ExN", "L6 IT ExN", "L6 CT ExN", "L6b ExN",
                                "Pvalb+ InN", "Sst+ InN", "Lamp5+ InN", "Vip+ InN",
                                "Micro", "Astro", "Oligo", "OPC", "Endo", "VLMC", "VSMC")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        text = element_text(family = "Arial")) +
  guides(color = guide_legend(title="Cell type"))


svglite("/xenium_Line61_final/ctx/ABC_comp/outs/ctx_ABC_comp_custom_genes.svg", 
        width = 8, height = 5.5)
plot(p2)
dev.off()
