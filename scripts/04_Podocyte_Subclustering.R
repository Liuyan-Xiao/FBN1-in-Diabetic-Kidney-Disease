
#pod subclustering

library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)

scRNA_ha <- readRDS("~/glomerulus_scdata/scRNA_ha_annotated.rds")

# 1. Extract podocytes
podo_cells <- subset(scRNA_ha, subset = celltype == "Podocyte")
table(podo_cells$orig.ident)
DefaultAssay(podo_cells) <- "RNA"
podo_cells <- CreateSeuratObject(
  counts = podo_cells@assays$RNA@counts,
  meta.data = podo_cells@meta.data
)

# 2. Normalize and identify variable features
podo_cells <- podo_cells %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData()

podo_cells_gene <- VariableFeatures(podo_cells)
anchor_features <- podo_cells_gene

# 3. PCA and Harmony
DefaultAssay(podo_cells) <- "RNA"
podo_cells <- RunPCA(
  object = podo_cells,
  npcs = 20,
  features = anchor_features
)

ElbowPlot(podo_cells)
table(podo_cells$orig.ident)

podo_cells <- RunHarmony(
  podo_cells,
  group.by.vars = "orig.ident",
  reduction = "pca",
  dims.use = 1:10
)

# 4. UMAP, t-SNE, and clustering
podo_cells <- RunUMAP(podo_cells, reduction = "harmony", dims = 1:10, seed.use = 20)
podo_cells <- RunTSNE(podo_cells, reduction = "harmony", dims = 1:10, seed.use = 20)
podo_cells <- FindNeighbors(podo_cells, reduction = "harmony", dims = 1:10)
podo_cells <- FindClusters(podo_cells, resolution = 1, verbose = FALSE)

Idents(podo_cells) <- "RNA_snn_res.1"
podo_cells$seurat_clusters <- podo_cells@active.ident
DefaultAssay(podo_cells) <- "RNA"

# 5. Identify cluster markers
markers_all <- FindAllMarkers(
  object = podo_cells,
  assay = "RNA",
  slot = "data",
  min.pct = 0.25,
  logfc.threshold = 0.25,
  only.pos = TRUE
)

# 6. Define podocyte-state gene signatures
injury_genes <- c("Des", "Clu", "C3", "Vcam1", "Tnfrsf12a",
  "Cxcl1", "B2m", "S100a6", "Lgals3")

stress_genes <- c("Fos", "Fosb", "Jun", "Junb", "Jund",
  "Egr1", "Atf3", "Dusp1", "Nr4a1", "Hspa1a", "Hspa1b")

mature_genes <- c("Nphs1", "Nphs2", "Synpo", "Podxl", "Ptpro", "Wt1")

# 7.Calculate module scores
podo_cells <- AddModuleScore(podo_cells, features = list(injury_genes), name = "InjuryScore")
podo_cells <- AddModuleScore(podo_cells, features = list(stress_genes), name = "StressScore")
podo_cells <- AddModuleScore(podo_cells, features = list(mature_genes), name = "MatureScore")

# 8. Annotate 
podo_cells_metadata <- podo_cells@meta.data
podo_cells_metadata$podo_state <- "Na"
cluster <- as.character(podo_cells_metadata$RNA_snn_res.1)

podo_cells_metadata$podo_state[cluster %in% c("0", "5", "6", "9")] <- "Stress_pod"
podo_cells_metadata$podo_state[cluster == "2"] <- "Unclass_pod"
podo_cells_metadata$podo_state[cluster %in% c("1", "8", "10")] <- "Injure_pod"
podo_cells_metadata$podo_state[cluster %in% c("3", "11", "7", "4")] <- "Mature_pod"

podo_cells@meta.data <- podo_cells_metadata

table(podo_cells$podo_state)

# 9. Visualize
VlnPlot(
  podo_cells,
  features = c("InjuryScore1", "StressScore1", "MatureScore1"),
  group.by = "podo_state"
  pt.size = 0,
  ncol = 3
)

DimPlot(
  podo_cells,
  reduction = "tsne",
  label = TRUE,
  split.by = "sample",
  group.by = "podo_state"
)


DotPlot(
  podo_cells,
  features = c(
    "Clu","Lgals3","Des","Spon2",
    "Fos", "Fosb", "Jun", "Egr1",
    "Nphs1", "Nphs2", "Synpo", "Wt1"),
  group.by = "podo_state") +
  scale_color_gradientn(colours = c("#FFFFFF","#FCE4EC","#F4B2CC",
                                    "#DE719A","#B83D72")
  )+RotatedAxis()



f <- FetchData(podo_cells, vars = "Fbn1")[,1]
podo_cells$Fbn1_status <- ifelse(
  f > 0,
  "Fbn1+",
  "Fbn1-"
)
DimPlot(
  podo_cells,
  reduction = "tsne",
  group.by = "Fbn1_status",
  split.by = "group",
  pt.size = 1,
  order = c( "Fbn1+","Fbn1-"),
  cols = c(
    "Fbn1-" = "grey90",
    "Fbn1+" = "#D73027"
  )
)
# 10. Save final object
saveRDS(podo_cells, "~/glomerulus_scdata/podo_cells.rds")


#FBN1 + cells(%)
df <- FetchData(
  podo_cells,
  vars = c("Fbn1", "podo_state")
)

df$Fbn1_pos <- df$Fbn1 > 0
table(df$podo_state,df$Fbn1_pos)

library(dplyr)
library(ggplot2)

state_pct <- df %>%
  group_by(podo_state) %>%
  summarise(
    pct_Fbn1 = mean(Fbn1_pos) * 100,
    .groups = "drop"
  )
#plot
state_pct <- df %>%
  group_by(podo_state) %>%
  summarise(
    Fbn1_pos_n = sum(Fbn1_pos),
    total_n = n(),
    pct_Fbn1 = Fbn1_pos_n / total_n * 100,
    .groups = "drop"
  ) %>%
  mutate(
    podo_state = factor(
      podo_state,
      levels = c("Injure_pod", "Stress_pod", "Mature_pod", "Unclass_pod")
    ),
    label = paste0(Fbn1_pos_n, "/", total_n, "\n", sprintf("%.2f%%", pct_Fbn1))
  )


ggplot(
  state_pct,
  aes(x = podo_state, y = pct_Fbn1, fill = podo_state)
) +
  geom_col(width = 0.5, color = "white", linewidth = 0.3) +
  geom_text(
    aes(label = label),
    vjust = -0.2,
    size = 4
  ) +
  scale_fill_manual(values = c(
    "Injure_pod"  = "#E99B68",
    "Stress_pod"  = "#8A91D4",
    "Mature_pod"  = "#79BFA5",
    "Unclass_pod" = "#D89AB6"
  )) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    x = NULL,
    y = "Fbn1+ cells (%)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black")
  )

