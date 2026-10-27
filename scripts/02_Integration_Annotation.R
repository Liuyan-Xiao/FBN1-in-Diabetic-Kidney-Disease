library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(patchwork)

setwd("~/glomerulus_scdata")
output_dir <- "~/glomerulus_scdata"

# 1. Load QC and DoubletFinder results
split_seurat <- readRDS("~/glomerulus_scdata/DoubletFinder_split_objects.rds")


# 2. Remove Doublets
for (i in seq_along(split_seurat)) {
  colnames(split_seurat[[i]]@meta.data)[ncol(split_seurat[[i]]@meta.data)] = "doublet_info"
  split_seurat[[i]] <- subset(split_seurat[[i]], subset = (doublet_info == "Singlet"))
}

# 3. SCTransform each sample
for (i in seq_along(split_seurat)) {
  split_seurat[[i]] <- SCTransform(
    split_seurat[[i]], assay = "RNA", seed.use = 1448145,
    vars.to.regress = "percent.mt", variable.features.n = 2000,
    return.only.var.genes = TRUE, verbose = TRUE
  )
}
saveRDS(split_seurat, file.path(output_dir, "Singlets_SCTransformed_split.rds"))

# 4. Merge samples
scRNA_ha <- merge(
  x = split_seurat[[1]],
  y = split_seurat[-1],
  project = "BTBR_glomerular_scRNAseq"
)

scRNA_ha$sample_id <- as.character(scRNA_ha$sample_id)
scRNA_ha$sample <- factor(scRNA_ha$sample, levels = sample_labels)
scRNA_ha$timepoint <- factor(scRNA_ha$timepoint, levels = c("control", "wk12", "wk21"))
scRNA_ha$group <- ifelse(scRNA_ha$sample_id %in% c("control_1", "control_2"),"ob/+","ob/ob")

# 5. PCA feature selection
anchor_features <- SelectIntegrationFeatures(split_seurat, nfeatures = 2000)


scRNA_ha <- RunPCA(
  scRNA_ha, npcs = 30, features = anchor_features,
  verbose = FALSE
)

# 6. Harmony batch correction
scRNA_ha <- RunHarmony(
  object = scRNA_ha,
  group.by.vars = "orig.ident",
  reduction = "pca",
  dims.use = 1:30,
  plot_convergence = TRUE
)

# 7. UMAP, t-SNE, and clustering
scRNA_ha <- RunUMAP(scRNA_ha, reduction = "harmony", dims = 1:30, seed.use = 42)
scRNA_ha <- RunTSNE(scRNA_ha, reduction = "harmony", dims = 1:30, seed.use = 42)
scRNA_ha <- FindNeighbors(scRNA_ha, reduction = "harmony", dims = 1:30, verbose = FALSE)
scRNA_ha <- FindClusters(scRNA_ha, resolution = seq(0.2, 1, by = 0.1), verbose = FALSE)


# 9. Find cluster markers
Idents(scRNA_ha) <- "SCT_snn_res.0.4"
scRNA_ha$seurat_clusters <- as.character(Idents(scRNA_ha))
DefaultAssay(scRNA_ha) <- "SCT"

scRNA_ha <- PrepSCTFindMarkers(scRNA_ha)
markers_all <- FindAllMarkers(
  scRNA_ha, min.pct = 0.25,
  logfc.threshold = 0.25, only.pos = TRUE
)

# 10. Cell-type annotation
metadata <- scRNA_ha@meta.data
metadata$celltype <- "Na"
cluster <- as.character(metadata$SCT_snn_res.0.4)

metadata$celltype[cluster %in% c("1","2","4","7","8","12")] <- "GEC"
metadata$celltype[cluster == "0"] <- "MC"
metadata$celltype[cluster == "3"] <- "Podocyte"
metadata$celltype[cluster %in% c("5","11")] <- "aEC"
metadata$celltype[cluster == "6"] <- "SMC"
metadata$celltype[cluster %in% c("9","10","13")] <- "Immune"
metadata$celltype[cluster == "14"] <- "PEC"
metadata$celltype[cluster %in% c("15","16","17","18","19")] <- "Tubular"

scRNA_ha@meta.data <- metadata
Idents(scRNA_ha) <- "celltype"

cell_colors <- c(
  "Tubular"     = "#A99ACB",
  "GEC"          = "#6BA8D7",
  "Immune" = "#DE9C6D",
  "MC"          = "#E47469",
  "PEC"         = "#F8D98B",
  "Podocyte"        = "#E6A0C5",
  "SMC"         = "#B2D990",
  "aEC"     = "#75B5A6"
)

DimPlot(
  scRNA_ha,
  reduction = "tsne",
  group.by = "celltype",raster = TRUE,
  raster.dpi = c(600, 600),
  cols = cell_colors,
  pt.size = 1.5,
  order = names(cell_colors)
) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    strip.text = element_text(
      size = 13,
      face = "bold"
    )
  )


# 11. Save final objects and metadata
saveRDS(scRNA_ha, file.path(output_dir, "scRNA_ha_annotated.rds"))
