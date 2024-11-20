library(Seurat)

xenium <- readRDS("/xenium_Line61_final/general/xenium_lognorm_integrated.rds")


# reading in csvs w/ cell IDs

nTg_1_cells <- read.csv("/xenium_Line61_final/xenium_upload1/ctx/nTg_1_cells.csv")
nTg_2_cells <- read.csv("/xenium_Line61_final/xenium_upload1/ctx/nTg_2_cells.csv")
nTg_3_cells <- read.csv("/xenium_Line61_final/xenium_upload1/ctx/nTg_3_cells.csv")
nTg_4_cells <- read.csv("/xenium_Line61_final/xenium_upload1/ctx/nTg_4_cells.csv")
Tg_1_cells <- read.csv("/xenium_Line61_final/xenium_upload1/ctx/Tg_1_cells.csv")
Tg_2_cells <- read.csv("/xenium_Line61_final/xenium_upload1/ctx/Tg_2_cells.csv")
Tg_3_cells <- read.csv("/xenium_Line61_final/xenium_upload1/ctx/Tg_3_cells.csv")
Tg_4_cells <- read.csv("/xenium_Line61_final/xenium_upload1/ctx/Tg_4_cells.csv")


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

xenium_ctx <- xenium[, cell_id_combined]

# re-run standard workflow

xenium_features <- rownames(xenium_ctx)

xenium_ctx <- ScaleData(xenium_ctx, features = xenium_features, vars.to.regress = c("nCount_Xenium"))
xenium_ctx <- RunPCA(xenium_ctx, features = xenium_features)
xenium_ctx <- RunUMAP(xenium_ctx, reduction = "harmony", dims = 1:30)
xenium_ctx <- FindNeighbors(xenium_ctx, reduction = "harmony", dims = 1:30)
xenium_ctx <- FindClusters(xenium_ctx, resolution = 1.2, cluster.name = "ctx.clusters")

xenium_ctx <- JoinLayers(xenium_ctx)

saveRDS(xenium_ctx, file = "/xenium_Line61_final/ctx/xenium_CTX-ONLY_lognorm_integrated_clustered.rds")
