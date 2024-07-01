library(SeuratDisk)
library(Seurat)
library(dplyr)

# convert adata to seurat using SeuratDisk
Convert("/new_10Xv2/raw_objects/WMB-10Xv2-Isocortex-1-raw.h5ad", "h5seurat")
Convert("/new_10Xv2/raw_objects/WMB-10Xv2-Isocortex-2-raw.h5ad", "h5seurat")
Convert("/new_10Xv2/raw_objects/WMB-10Xv2-Isocortex-3-raw.h5ad", "h5seurat")
Convert("/new_10Xv2/raw_objects/WMB-10Xv2-Isocortex-4-raw.h5ad", "h5seurat")


# load new h5seurat objects
isoctx_1 <- LoadH5Seurat("/new_10Xv2/raw_objects/WMB-10Xv2-Isocortex-1-raw.h5seurat", meta.data = F)
isoctx_2 <- LoadH5Seurat("/new_10Xv2/raw_objects/WMB-10Xv2-Isocortex-2-raw.h5seurat", meta.data = F)
isoctx_3 <- LoadH5Seurat("/new_10Xv2/raw_objects/WMB-10Xv2-Isocortex-3-raw.h5seurat", meta.data = F)
isoctx_4 <- LoadH5Seurat("/new_10Xv2/raw_objects/WMB-10Xv2-Isocortex-4-raw.h5seurat", meta.data = F)


# load cell-level metadata from ABC
metadata <- read.csv("/metadata/cell_metadata_with_cluster_annotation.csv")
meta.df <- as.data.frame(metadata)


# get cell IDs from seurat object to subset metadata for only the cells contained in the seurat object
cells1 <- Cells(isoctx_1)
meta.df1 <- meta.df[meta.df$cell_label %in% cells1, ]

cells2 <- Cells(isoctx_2)
meta.df2 <- meta.df[meta.df$cell_label %in% cells2, ]

cells3 <- Cells(isoctx_3)
meta.df3 <- meta.df[meta.df$cell_label %in% cells3, ]

cells4 <- Cells(isoctx_4)
meta.df4 <- meta.df[meta.df$cell_label %in% cells4, ]

# add metadata to seurat objects

isoctx_1[['region_of_interest_acronym']] <- meta.df1$region_of_interest_acronym[match(rownames(isoctx_1@meta.data), meta.df1$cell_label)]
isoctx_1[['donor_label']] <- meta.df1$donor_label[match(rownames(isoctx_1@meta.data), meta.df1$cell_label)]
isoctx_1[['donor_sex']] <- meta.df1$donor_sex[match(rownames(isoctx_1@meta.data), meta.df1$cell_label)]
isoctx_1[['neurotransmitter']] <- meta.df1$neurotransmitter[match(rownames(isoctx_1@meta.data), meta.df1$cell_label)]
isoctx_1[['division']] <- meta.df1$division[match(rownames(isoctx_1@meta.data), meta.df1$cell_label)]
isoctx_1[['class']] <- meta.df1$class[match(rownames(isoctx_1@meta.data), meta.df1$cell_label)]
isoctx_1[['subclass']] <- meta.df1$subclass[match(rownames(isoctx_1@meta.data), meta.df1$cell_label)]

isoctx_2[['region_of_interest_acronym']] <- meta.df2$region_of_interest_acronym[match(rownames(isoctx_2@meta.data), meta.df2$cell_label)]
isoctx_2[['donor_label']] <- meta.df2$donor_label[match(rownames(isoctx_2@meta.data), meta.df2$cell_label)]
isoctx_2[['donor_sex']] <- meta.df2$donor_sex[match(rownames(isoctx_2@meta.data), meta.df2$cell_label)]
isoctx_2[['neurotransmitter']] <- meta.df2$neurotransmitter[match(rownames(isoctx_2@meta.data), meta.df2$cell_label)]
isoctx_2[['division']] <- meta.df2$division[match(rownames(isoctx_2@meta.data), meta.df2$cell_label)]
isoctx_2[['class']] <- meta.df2$class[match(rownames(isoctx_2@meta.data), meta.df2$cell_label)]
isoctx_2[['subclass']] <- meta.df2$subclass[match(rownames(isoctx_2@meta.data), meta.df2$cell_label)]

isoctx_3[['region_of_interest_acronym']] <- meta.df3$region_of_interest_acronym[match(rownames(isoctx_3@meta.data), meta.df3$cell_label)]
isoctx_3[['donor_label']] <- meta.df3$donor_label[match(rownames(isoctx_3@meta.data), meta.df3$cell_label)]
isoctx_3[['donor_sex']] <- meta.df3$donor_sex[match(rownames(isoctx_3@meta.data), meta.df3$cell_label)]
isoctx_3[['neurotransmitter']] <- meta.df3$neurotransmitter[match(rownames(isoctx_3@meta.data), meta.df3$cell_label)]
isoctx_3[['division']] <- meta.df3$division[match(rownames(isoctx_3@meta.data), meta.df3$cell_label)]
isoctx_3[['class']] <- meta.df3$class[match(rownames(isoctx_3@meta.data), meta.df3$cell_label)]
isoctx_3[['subclass']] <- meta.df3$subclass[match(rownames(isoctx_3@meta.data), meta.df3$cell_label)]

isoctx_4[['region_of_interest_acronym']] <- meta.df4$region_of_interest_acronym[match(rownames(isoctx_4@meta.data), meta.df4$cell_label)]
isoctx_4[['donor_label']] <- meta.df4$donor_label[match(rownames(isoctx_4@meta.data), meta.df4$cell_label)]
isoctx_4[['donor_sex']] <- meta.df4$donor_sex[match(rownames(isoctx_4@meta.data), meta.df4$cell_label)]
isoctx_4[['neurotransmitter']] <- meta.df4$neurotransmitter[match(rownames(isoctx_4@meta.data), meta.df4$cell_label)]
isoctx_4[['division']] <- meta.df4$division[match(rownames(isoctx_4@meta.data), meta.df4$cell_label)]
isoctx_4[['class']] <- meta.df4$class[match(rownames(isoctx_4@meta.data), meta.df4$cell_label)]
isoctx_4[['subclass']] <- meta.df4$subclass[match(rownames(isoctx_4@meta.data), meta.df4$cell_label)]


# pull cell counts by division and subclass

Idents(isoctx_1) <- "subclass"
subclass_cts_by_sample <- table(Idents(isoctx_1), isoctx_1$donor_label)
write.csv(subclass_cts_by_sample,
          file = "/new_10Xv2/outs/isoctx_1_subclass_ctsbysample.csv")

Idents(isoctx_2) <- "subclass"
subclass_cts_by_sample2 <- table(Idents(isoctx_2), isoctx_2$donor_label)
write.csv(subclass_cts_by_sample2,
          file = "/new_10Xv2/outs/isoctx_2_subclass_ctsbysample.csv")

Idents(isoctx_3) <- "subclass"
subclass_cts_by_sample3 <- table(Idents(isoctx_3), isoctx_3$donor_label)
write.csv(subclass_cts_by_sample3,
          file = "/new_10Xv2/outs/isoctx_3_subclass_ctsbysample.csv")

Idents(isoctx_4) <- "subclass"
subclass_cts_by_sample4 <- table(Idents(isoctx_4), isoctx_4$donor_label)
write.csv(subclass_cts_by_sample4,
          file = "/new_10Xv2/outs/isoctx_4_subclass_ctsbysample.csv")


# filter for desired animals
animals_keep <- c("Snap25-IRES2-Cre;Ai14-352353",
                  "Snap25-IRES2-Cre;Ai14-352356",
                  "Snap25-IRES2-Cre;Ai14-352357",
                  "Snap25-IRES2-Cre;Ai14-355878",
                  "Snap25-IRES2-Cre;Ai14-365617",
                  "Snap25-IRES2-Cre;Ai14-365619",
                  "Snap25-IRES2-Cre;Ai14-366676",
                  "Snap25-IRES2-Cre;Ai14-366678",
                  "Snap25-IRES2-Cre;Ai14-371230",
                  "Snap25-IRES2-Cre;Ai14-371232",
                  "Snap25-IRES2-Cre;Ai14-372312",
                  "Snap25-IRES2-Cre;Ai14-372317",
                  "Snap25-IRES2-Cre;Ai14-373822",
                  "Snap25-IRES2-Cre;Ai14-374158",
                  "Snap25-IRES2-Cre;Ai14-374160",
                  "Snap25-IRES2-Cre;Ai14-381296",
                  "Snap25-IRES2-Cre;Ai14-395345",
                  "Snap25-IRES2-Cre;Ai14-374160",
                  "Snap25-IRES2-Cre;Ai14-410107",
                  "Snap25-IRES2-Cre;Ai14-410108",
                  "Snap25-IRES2-Cre;Ai14-412204",
                  "Snap25-IRES2-Cre;Ai14-412866",
                  "Snap25-IRES2-Cre;Ai14-414902",
                  "Snap25-IRES2-Cre;Ai14-414905",
                  "Snap25-IRES2-Cre;Ai14-414906",
                  "Snap25-IRES2-Cre;Ai14-446716")

Idents(isoctx_1) <- "donor_label"
isoctx_1 <- subset(isoctx_1, idents = animals_keep)
Idents(isoctx_2) <- "donor_label"
isoctx_2 <- subset(isoctx_2, idents = animals_keep)
Idents(isoctx_3) <- "donor_label"
isoctx_3 <- subset(isoctx_3, idents = animals_keep)
Idents(isoctx_4) <- "donor_label"
isoctx_4 <- subset(isoctx_4, idents = animals_keep)


# downsample each set

set.seed(2)
isoctx_1_downsampled <- isoctx_1[, sample(colnames(isoctx_1), size = 60000, replace=F)]
isoctx_2_downsampled <- isoctx_2[, sample(colnames(isoctx_2), size = 60000, replace=F)]
isoctx_3_downsampled <- isoctx_3[, sample(colnames(isoctx_3), size = 60000, replace=F)]
isoctx_4_downsampled <- isoctx_4[, sample(colnames(isoctx_4), size = 60000, replace=F)]


# merge downsampled objects

isoctx <- merge(x = isoctx_1_downsampled, y = c(isoctx_2_downsampled, isoctx_3_downsampled, isoctx_4_downsampled))


# celltype counts on whole object

Idents(isoctx) <- "subclass"
subclass_cts_by_sample_all <- table(Idents(isoctx), isoctx$donor_label)
write.csv(subclass_cts_by_sample_all,
          file = "/new_10Xv2/outs/isoctx_all_subclass_ctsbysample.csv")



# save raw merged file

saveRDS(isoctx, file = "/new_10Xv2/isoctx_merged_240k_raw.rds")

isoctx <- readRDS("/new_10Xv2/isoctx_merged_240k_raw.rds")


# filtering further 

Idents(isoctx) <- "division"
isoctx <- subset(isoctx, idents = c("1 Pallium glutamatergic",
                                    "2 Subpallium GABAergic",
                                    "5 Neuroglial",
                                    "6 Vascular",
                                    "7 Immune"))

# standard seurat workflow for visualizing clusters

isoctx <- NormalizeData(isoctx, normalization.method = "LogNormalize", scale.factor = 1e6, assay = "RNA")
isoctx <- FindVariableFeatures(isoctx, selection.method = "vst", nfeatures = 2000)
isoctx <- ScaleData(object = isoctx)
isoctx <- RunPCA(object = isoctx)
isoctx <- FindNeighbors(object = isoctx, dims = 1:30)
isoctx <- RunUMAP(object = isoctx, dims = 1:30)


# adding "celltype" column in metadata

Idents(isoctx) <- "subclass"
levels(isoctx)

isoctx$celltype <- ifelse(isoctx$subclass == "007 L2/3 IT CTX Glut","L2/3 IT ExN",
                          ifelse(isoctx$subclass == "006 L4/5 IT CTX Glut", "L4/5 IT ExN",
                                 ifelse(isoctx$subclass == "005 L5 IT CTX Glut", "L5 IT ExN",
                                        ifelse(isoctx$subclass == "022 L5 ET CTX Glut", "L5 ET ExN",
                                               ifelse(isoctx$subclass == "030 L5 NP CTX Glut", "L5 NP ExN",
                                                      ifelse(isoctx$subclass == "004 L6 IT CTX Glut", "L6 IT ExN",
                                                             ifelse(isoctx$subclass == "028 L6 CT CTX Glut", "L6 CT ExN",
                                                                    ifelse(isoctx$subclass == "027 L6b CTX Glut", "L6b ExN",
                                                                           ifelse(isoctx$subclass == "042 Pvalb Gaba", "Pvalb+ InN",
                                                                                  ifelse(isoctx$subclass == "036 Vip Gaba", "Vip+ InN",
                                                                                         ifelse(isoctx$subclass == "039 Lamp5 Gaba", "Lamp5+ InN",
                                                                                                ifelse(isoctx$subclass == "043 Sst Gaba", "Sst+ InN",
                                                                                                       ifelse(isoctx$subclass == "302 Microglia NN", "Micro",
                                                                                                              ifelse(isoctx$subclass == "287 Astro-TE NN", "Astro",
                                                                                                                     ifelse(isoctx$subclass == "295 Oligo NN", "Oligo",
                                                                                                                            ifelse(isoctx$subclass == "294 OPC NN", "OPC",
                                                                                                                                   ifelse(isoctx$subclass == "299 Peri NN", "Peri",
                                                                                                                                          ifelse(isoctx$subclass == "300 SMC NN", "VSMC",
                                                                                                                                                 ifelse(isoctx$subclass == "298 VLMC NN", "VLMC",
                                                                                                                                                        ifelse(isoctx$subclass == "301 Endo NN", "Endo", "Other"))))))))))))))))))))


Idents(isoctx) <- "celltype"
levels(isoctx)
levels(isoctx) <- c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", "L5 ET ExN", "L5 NP ExN", "L6 IT ExN", "L6 CT ExN", "L6b ExN",
                    "Pvalb+ InN", "Sst+ InN", "Vip+ InN", "Lamp5+ InN", "Micro", "Astro", "Oligo", "OPC", "Peri", 
                    "VSMC", "VLMC", "Endo", "Other")


saveRDS(isoctx, "/new_10Xv2/isocortex_all_lognorm_filtered_seurat.rds")
