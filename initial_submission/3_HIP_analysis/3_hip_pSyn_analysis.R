library(Seurat)

xenium_hip <- readRDS("/xenium_Line61_final/hip/xenium_HIP_metadata.rds")


# overlap analysis w/ polygons from QuPath
library(sf)
library(dplyr)
library(tidyverse)
library(ggplot2)

### This is everything written out for one sample, for the other 3 it will just be the code

### Tg_1 ###
# read in centroids of pSyn polygons from QuPath & turn into sf object

tg_1_pSyn_centroids <- read.csv("/xenium_Line61_final/xenium_upload1/general/tg_1_centroids.csv")
tg_1_pSyn.df <- as.data.frame(tg_1_pSyn_centroids)
tg_1_pSyn.sf <- st_as_sf(tg_1_pSyn.df, coords = c("Centroid.X.µm", "Centroid.Y.µm"))

# extract polygons from xenium cells and convert to sf object
# polygons for xenium cells are stored at xenium@images$fov@boundaries$segmentation@polygons

tg_1_cell_polygons <- xenium_hip@images$fov.5@boundaries$segmentation@polygons
# extract cell names
tg_1_cell_names <- names(tg_1_cell_polygons)
# extract coords for individual cells
tg_1_cell_coords <- lapply(tg_1_cell_names, function(cell_name) {
  coords <- tg_1_cell_polygons[[cell_name]]@Polygons[[1]]@coords
  data.frame(cell_name = cell_name, x = coords[,1], y = coords[,2])
})
# combine cell coords to single long-form dataframe
tg_1_cell_coords.df <- do.call(rbind, tg_1_cell_coords)

# function for converting cell coords df to list of polygons in sfg format
convert_to_sf <- function(df, id_column, x_column, y_column, crs = 4326) {
  polygons <- list()
  
  for (id in unique(df[[id_column]])) {
    cell_df <- df[df[[id_column]] == id, ]
    coords <- cbind(cell_df[[x_column]], cell_df[[y_column]])
    poly <- st_polygon(list(coords))
    polygons[[id]] <- poly
  }
  
  return(polygons)
}

# run the function
tg_1_polygons.sf <- convert_to_sf(df = tg_1_cell_coords.df, id_column = "cell_name", x_column = "x", y_column = "y")

# convert list of sfg polygons to big sfc object
tg_1_polygons.sfc <- st_as_sfc(tg_1_polygons.sf)

# use st_contains to retrieve cells w/ pSyn
tg_1_pSyn_cells <- st_contains(tg_1_polygons.sfc, tg_1_pSyn.sf)

# make dataframe w/ cell IDs & whether or not they have pSyn
tg_1_pSyn_cells.df <- as.data.frame(tg_1_pSyn_cells)
tg_1_cells_w_pSyn.df <- as.data.frame(names(tg_1_polygons.sfc))
tg_1_cells_w_pSyn.df$has_pSyn <- rownames(tg_1_cells_w_pSyn.df) %in% tg_1_pSyn_cells.df$row.id
colnames(tg_1_cells_w_pSyn.df) <- c("cell_ID", "has_pSyn")


# can plot the 2 sf objects to confirm overlap
ggplot() +
  geom_sf(data = tg_1_polygons.sfc, fill = "transparent", color = "blue") +
  geom_sf(data = tg_1_pSyn.sf, color = "red", size = 0.5) +
  theme_minimal()




### Tg_2 ###

tg_2_pSyn_centroids <- read.csv("/xenium_Line61_final/xenium_upload1/general/tg_2_centroids.csv")
tg_2_pSyn.df <- as.data.frame(tg_2_pSyn_centroids)
tg_2_pSyn.sf <- st_as_sf(tg_2_pSyn.df, coords = c("Centroid.X.µm", "Centroid.Y.µm"))
tg_2_cell_polygons <- xenium_hip@images$fov.6@boundaries$segmentation@polygons
tg_2_cell_names <- names(tg_2_cell_polygons)
tg_2_cell_coords <- lapply(tg_2_cell_names, function(cell_name) {
  coords <- tg_2_cell_polygons[[cell_name]]@Polygons[[1]]@coords
  data.frame(cell_name = cell_name, x = coords[,1], y = coords[,2])
})
tg_2_cell_coords.df <- do.call(rbind, tg_2_cell_coords)
tg_2_polygons.sf <- convert_to_sf(df = tg_2_cell_coords.df, id_column = "cell_name", x_column = "x", y_column = "y")
tg_2_polygons.sfc <- st_as_sfc(tg_2_polygons.sf)
tg_2_pSyn_cells <- st_contains(tg_2_polygons.sfc, tg_2_pSyn.sf)
tg_2_pSyn_cells.df <- as.data.frame(tg_2_pSyn_cells)
tg_2_cells_w_pSyn.df <- as.data.frame(names(tg_2_polygons.sfc))
tg_2_cells_w_pSyn.df$has_pSyn <- rownames(tg_2_cells_w_pSyn.df) %in% tg_2_pSyn_cells.df$row.id
colnames(tg_2_cells_w_pSyn.df) <- c("cell_ID", "has_pSyn")



### Tg_3 ###

tg_3_pSyn_centroids <- read.csv("/xenium_Line61_final/xenium_upload1/general/tg_3_centroids.csv")
tg_3_pSyn.df <- as.data.frame(tg_3_pSyn_centroids)
tg_3_pSyn.sf <- st_as_sf(tg_3_pSyn.df, coords = c("Centroid.X.µm", "Centroid.Y.µm"))
tg_3_cell_polygons <- xenium_hip@images$fov.7@boundaries$segmentation@polygons
tg_3_cell_names <- names(tg_3_cell_polygons)
tg_3_cell_coords <- lapply(tg_3_cell_names, function(cell_name) {
  coords <- tg_3_cell_polygons[[cell_name]]@Polygons[[1]]@coords
  data.frame(cell_name = cell_name, x = coords[,1], y = coords[,2])
})
tg_3_cell_coords.df <- do.call(rbind, tg_3_cell_coords)
tg_3_polygons.sf <- convert_to_sf(df = tg_3_cell_coords.df, id_column = "cell_name", x_column = "x", y_column = "y")
tg_3_polygons.sfc <- st_as_sfc(tg_3_polygons.sf)
tg_3_pSyn_cells <- st_contains(tg_3_polygons.sfc, tg_3_pSyn.sf)
tg_3_pSyn_cells.df <- as.data.frame(tg_3_pSyn_cells)
tg_3_cells_w_pSyn.df <- as.data.frame(names(tg_3_polygons.sfc))
tg_3_cells_w_pSyn.df$has_pSyn <- rownames(tg_3_cells_w_pSyn.df) %in% tg_3_pSyn_cells.df$row.id
colnames(tg_3_cells_w_pSyn.df) <- c("cell_ID", "has_pSyn")



### Tg_4 ###

tg_4_pSyn_centroids <- read.csv("/xenium_Line61_final/xenium_upload1/general/tg_4_centroids.csv")
tg_4_pSyn.df <- as.data.frame(tg_4_pSyn_centroids)
tg_4_pSyn.sf <- st_as_sf(tg_4_pSyn.df, coords = c("Centroid.X.µm", "Centroid.Y.µm"))
tg_4_cell_polygons <- xenium_hip@images$fov.8@boundaries$segmentation@polygons
tg_4_cell_names <- names(tg_4_cell_polygons)
tg_4_cell_coords <- lapply(tg_4_cell_names, function(cell_name) {
  coords <- tg_4_cell_polygons[[cell_name]]@Polygons[[1]]@coords
  data.frame(cell_name = cell_name, x = coords[,1], y = coords[,2])
})
tg_4_cell_coords.df <- do.call(rbind, tg_4_cell_coords)
tg_4_polygons.sf <- convert_to_sf(df = tg_4_cell_coords.df, id_column = "cell_name", x_column = "x", y_column = "y")
tg_4_polygons.sfc <- st_as_sfc(tg_4_polygons.sf)
tg_4_pSyn_cells <- st_contains(tg_4_polygons.sfc, tg_4_pSyn.sf)
tg_4_pSyn_cells.df <- as.data.frame(tg_4_pSyn_cells)
tg_4_cells_w_pSyn.df <- as.data.frame(names(tg_4_polygons.sfc))
tg_4_cells_w_pSyn.df$has_pSyn <- rownames(tg_4_cells_w_pSyn.df) %in% tg_4_pSyn_cells.df$row.id
colnames(tg_4_cells_w_pSyn.df) <- c("cell_ID", "has_pSyn")


# merge dataframes w/ cells w/ pSyn

all_cells_w_pSyn.df <- rbind(tg_1_cells_w_pSyn.df, tg_2_cells_w_pSyn.df, tg_3_cells_w_pSyn.df, tg_4_cells_w_pSyn.df)


# add in pSyn status to metadata of Seurat object

all_cells_w_pSyn.df <- column_to_rownames(all_cells_w_pSyn.df, "cell_ID")
matching_cells <- intersect(colnames(xenium_hip), rownames(all_cells_w_pSyn.df))
xenium_hip$has_pSyn <- NA
xenium_hip$has_pSyn[matching_cells] <- all_cells_w_pSyn.df[matching_cells, "has_pSyn"]


# save Xenium object

saveRDS(xenium_hip, file = "/xenium_Line61_final/hip/xenium_HIP_metadata_pSyn.rds")
