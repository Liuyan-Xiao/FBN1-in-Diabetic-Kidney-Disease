library(readr)
library(dplyr)
library(tibble)
library(stringr)
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(patchwork)

# read file
file_path <- "~/data/GSE307956_core_table_gene.txt.gz"
df <- read_tsv(file_path)

Counts <- df[,c(3,22:30)]
colnames(Counts)[1] <- "gene_symbol"
Counts <- na.omit(Counts)


Counts <- column_to_rownames(Counts, var = "gene_symbol")

count_matrix <- as.matrix(Counts)
count_matrix <- round(as.matrix(count_matrix))
colnames(count_matrix)

# 1.sample name
sample_names <- colnames(count_matrix)

# 2.group
group <- case_when(
  str_detect(sample_names, "read_count_N") ~ "Normal",
  str_detect(sample_names, "read_count_M") ~ "Middle",
  str_detect(sample_names, "read_count_H") ~ "High"
)

# 3. 
colData <- data.frame(
  row.names = sample_names,
  group = factor(group,
                 levels = c("Normal", "Middle", "High")))


dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = colData,
  design = ~ group
)
# filter
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# 
dds <- DESeq(dds)


res_high <- results(dds, contrast = c("group", "High", "Normal"))
res_middle <- results(dds, contrast = c("group", "Middle", "Normal"))

# PCA using the variance-stabilized expression matrix.
vsd <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
colors <- c(Normal = "#6BA8D7", Middle = "#F8D98B", High = "#E47469")

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = group, label = name)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(vjust = -1.2, size = 3) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  scale_color_manual(values = colors) +
  theme_classic(base_size = 14) +
  theme(legend.title = element_blank())

# FBN1 abundance: log2(DESeq2-normalized counts + 1).
norm_counts <- counts(dds, normalized = TRUE)
plot_df <- data.frame(
  sample = colnames(norm_counts),
  group = as.character(colData(dds)$group),
  expression = as.numeric(log2(norm_counts["FBN1", ] + 1))
)


stat_df <- data.frame(
  group1 = c("Normal", "Normal"),
  group2 = c("Middle", "High"),
  FDR = c(res_middle["FBN1", "padj"], res_high["FBN1", "padj"])
)
stat_df$label <- ifelse(
  is.na(stat_df$FDR), "FDR = NA",
  paste0("FDR = ", formatC(stat_df$FDR, format = "f", digits = 3))
)
ymax <- max(plot_df$expression)
stat_df$y.position <- c(ymax + 0.3, ymax + 0.8)

p_fbn1 <- ggplot(plot_df, aes(group, expression, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.8,
               color = "#333333", linewidth = 0.6) +
  geom_jitter(width = 0.08, size = 2.8, color = "#222222") +
  stat_pvalue_manual(stat_df, label = "label", xmin = "group1", xmax = "group2",
                     y.position = "y.position", tip.length = 0.01,
                     bracket.size = 0.5, size = 4, inherit.aes = FALSE) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.25))) +
  labs(title = "FBN1 Expression", x = NULL,
       y = "FBN1 expression (log2 normalized counts + 1)") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text = element_text(color = "black"),
        legend.position = "none")

p_pca+p_fbn1



