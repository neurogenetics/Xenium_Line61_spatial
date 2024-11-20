library(Seurat)

xenium <- readRDS("/xenium_Line61_final/general/xenium_lognorm_integrated.rds")


# reading in csvs w/ cell IDs

nTg_1_cells <- read.csv("/xenium_Line61_final/xenium_upload1/hip/nTg_1_HIP_cells.csv")
nTg_2_cells <- read.csv("/xenium_Line61_final/xenium_upload1/hip/nTg_2_HIP_cells.csv")
nTg_3_cells <- read.csv("/xenium_Line61_final/xenium_upload1/hip/nTg_3_HIP_cells.csv")
nTg_4_cells <- read.csv("/xenium_Line61_final/xenium_upload1/hip/nTg_4_HIP_cells.csv")
Tg_1_cells <- read.csv("/xenium_Line61_final/xenium_upload1/hip/Tg_1_HIP_cells.csv")
Tg_2_cells <- read.csv("/xenium_Line61_final/xenium_upload1/hip/Tg_2_HIP_cells.csv")
Tg_3_cells <- read.csv("/xenium_Line61_final/xenium_upload1/hip/Tg_3_HIP_cells.csv")
Tg_4_cells <- read.csv("/xenium_Line61_final/xenium_upload1/hip/Tg_4_HIP_cells.csv")


# add _x to cell IDs for whichever sample 

nTg_1_cells$Cell.ID <- paste0(nTg_1_cells$Cell.ID, "_1")
nTg_2_cells$Cell.ID <- paste0(nTg_2_cells$Cell.ID, "_2")
nTg_3_cells$Cell.ID <- paste0(nTg_3_cells$Cell.ID, "_3")
nTg_4_cells$Cell.ID <- paste0(nTg_4_cells$Cell.ID, "_4")
Tg_1_cells$Cell.ID <- paste0(Tg_1_cells$Cell.ID, "_5")
Tg_2_cells$Cell.ID <- paste0(Tg_2_cells$Cell.ID, "_6")
Tg_3_cells$Cell.ID <- paste0(Tg_3_cells$Cell.ID, "_7")
Tg_4_cells$Cell.ID <- paste0(Tg_4_cells$Cell.ID, "_8")


nTg_1_cells.df <- as.data.frame(nTg_1_cells)
nTg_2_cells.df <- as.data.frame(nTg_2_cells)
nTg_3_cells.df <- as.data.frame(nTg_3_cells)
nTg_4_cells.df <- as.data.frame(nTg_4_cells)
Tg_1_cells.df <- as.data.frame(Tg_1_cells)
Tg_2_cells.df <- as.data.frame(Tg_2_cells)
Tg_3_cells.df <- as.data.frame(Tg_3_cells)
Tg_4_cells.df <- as.data.frame(Tg_4_cells)

df_list <- list(nTg_1_cells.df, nTg_2_cells.df, nTg_3_cells.df, nTg_4_cells.df,
                Tg_1_cells.df, Tg_2_cells.df, Tg_3_cells.df, Tg_4_cells.df)

cell_id_list <- lapply(df_list, function(df) df$Cell.ID)

# list of all cortical cell IDs

cell_id_combined <- unlist(cell_id_list)


# subset Seurat obj based on list of Cell ID

xenium_hip <- xenium[, cell_id_combined]

# re-run standard workflow

xenium_features <- rownames(xenium_hip)

xenium_hip <- ScaleData(xenium_hip, features = xenium_features, vars.to.regress = c("nCount_Xenium"))
xenium_hip <- RunPCA(xenium_hip, features = xenium_features)
xenium_hip <- RunUMAP(xenium_hip, reduction = "harmony", dims = 1:30)
xenium_hip <- FindNeighbors(xenium_hip, reduction = "harmony", dims = 1:30)
xenium_hip <- FindClusters(xenium_hip, resolution = 0.3, cluster.name = "hip.clusters")

xenium_hip <- JoinLayers(xenium_hip)

saveRDS(xenium_hip, file = "/xenium_Line61_final/hip/xenium_HIP-ONLY_lognorm_integrated_clustered.rds")
