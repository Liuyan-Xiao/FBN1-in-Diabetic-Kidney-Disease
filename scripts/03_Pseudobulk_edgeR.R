## Pseudobulk analysis of Fbn1 in podocytes
library(Seurat)
library(edgeR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# 1. Aggregate counts by biological sample
scRNA_ha <- readRDS("~/glomerulus_scdata/scRNA_ha_annotated.rds")
podo_cells <- subset(scRNA_ha, subset = celltype == "Podocyte")
table(podo_cells$orig.ident)
DefaultAssay(podo_cells) <- "RNA"
podo_cells <- CreateSeuratObject(
  counts = podo_cells@assays$RNA@counts,
  meta.data = podo_cells@meta.data
)

pb <- AggregateExpression(
  podo_cells,
  assays = "RNA",
  group.by = "sample",
  slot = "counts",
  return.seurat = FALSE
)$RNA

# 2. Prepare sample metadata
meta_pb <- podo_cells@meta.data %>%
  distinct(sample, timepoint) %>%
  filter(sample %in% colnames(pb))

meta_pb <- meta_pb[match(colnames(pb), meta_pb$sample), ]


meta_pb$condition <- factor(
  ifelse(meta_pb$timepoint == "control", "ob/+", "ob/ob"),
  levels = c("ob/+", "ob/ob")
)

# 3. edgeR differential expression analysis
dge <- DGEList(counts = pb,group = meta_pb$condition)
keep <- filterByExpr(dge, group = meta_pb$condition)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

design <- model.matrix(~ condition, data = meta_pb)

dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

qlf <- glmQLFTest(
  fit,
  coef = "conditionob/ob"
)
res <- topTags(qlf, n = Inf)$table
# Extract Fbn1 differential expression results
fbn1_res <- res["Fbn1", ]
print(fbn1_res)

# 4. Extract sample-level Fbn1 expression
dge_norm <- calcNormFactors(DGEList(counts = pb))

logCPM <- cpm(dge_norm,log = TRUE,prior.count = 1)

fbn1_pb <- data.frame(
  sample = colnames(logCPM),
  Fbn1 = as.numeric(logCPM["Fbn1", ])
) %>%
  left_join(meta_pb, by = "sample")

fbn1_pb$condition <- factor(fbn1_pb$condition,levels = c("ob/+", "ob/ob"))

# 5. Statistical annotation
fdr_value <- fbn1_res$FDR
stat_df <- data.frame(
  group1 = "ob/+",
  group2 = "ob/ob",
  y.position = max(fbn1_pb$Fbn1, na.rm = TRUE) + 0.3,
  label = paste0("FDR = ", formatC(fdr_value, format = "f", digits = 3))
)

# 6. Plot sample-level pseudobulk expression
ggplot(
  fbn1_pb,
  aes(x = condition, y = Fbn1, fill = condition)
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_point(
    size = 2,
    shape = 21,
    position = position_jitter(width = 0.07)
  ) +
  stat_pvalue_manual(
    stat_df,
    label = "label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "ob/+" = "#8CCDB7",
      "ob/ob" = "#bb5256"
    )
  ) +
  labs(
    x = NULL,
    y = "Fbn1 pseudobulk expression (log2 CPM)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10)
  )
