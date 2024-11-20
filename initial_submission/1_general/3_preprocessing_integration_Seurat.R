library(Seurat)

# load in Xenium objects

nTg_1 <- LoadXenium("/xenium_5um_expansion/nTg_1/outs", fov = "fov")
nTg_2 <- LoadXenium("/xenium_5um_expansion/nTg_2/outs", fov = "fov")
nTg_3 <- LoadXenium("/xenium_5um_expansion/nTg_3/outs", fov = "fov")
nTg_4 <- LoadXenium("/xenium_5um_expansion/nTg_4/outs", fov = "fov")
Tg_1 <- LoadXenium("/xenium_5um_expansion/Tg_1/outs", fov = "fov")
Tg_2 <- LoadXenium("/xenium_5um_expansion/Tg_2/outs", fov = "fov")
Tg_3 <- LoadXenium("/xenium_5um_expansion/Tg_3/outs", fov = "fov")
Tg_4 <- LoadXenium("/xenium_5um_expansion/Tg_4/outs", fov = "fov")


# merge objects & split by orig.ident

xenium <- merge(x = nTg_1, y = c(nTg_2, nTg_3, nTg_4, Tg_1, Tg_2, Tg_3, Tg_4))
xenium <- subset(xenium, nCount_Xenium > 0)


# normalization / integration

xenium <- NormalizeData(xenium)
# don't need to run FindVariableFeatures, just use all features for this
xenium <- ScaleData(xenium, features = rownames(xenium), vars.to.regress = c("nCount_Xenium"))
xenium.features <- rownames(xenium)
xenium <- RunPCA(xenium, features = xenium.features)

xenium <- IntegrateLayers(
  object = xenium, method = HarmonyIntegration, features = xenium.features,
  orig.reduction = "pca", new.reduction = "harmony", verbose = T
)

saveRDS(xenium, file = "/xenium_Line61_final/general/xenium_lognorm_integrated.rds")
