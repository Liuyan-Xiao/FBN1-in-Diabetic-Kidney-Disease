library(Seurat)
library(data.table)
library(Matrix)
library(DoubletFinder)
library(ggplot2)
library(patchwork)

setwd("~/glomerulus_scdata")
output_dir <- "~/glomerulus_scdata"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

files <- c("GSM4409518_BTBR_control_1.txt.gz","GSM4409519_BTBR_control_2.txt.gz","GSM4409520_BTBR_wk12_1.txt.gz","GSM4409521_BTBR_wk12_2.txt.gz","GSM4409522_BTBR_wk21_1.txt.gz","GSM4409523_BTBR_wk21_2.txt.gz")
sample_ids <- c("control_1","control_2","wk12_1","wk12_2","wk21_1","wk21_2")
sample_labels <- c("ob/+_1","ob/+_2","ob/ob_1","ob/ob_2","ob/ob_3","ob/ob_4")
input_paths <- file.path("~/glomerulus_scdata", files)

# 1. Read data and create Seurat objects
seurat_objects <- list()
for(i in seq_along(files)){
  counts <- fread(input_paths[i],data.table=FALSE)
  gene_names <- as.character(counts$SYMBOL)
  idx <- is.na(gene_names) | gene_names==""
  gene_names[idx] <- as.character(counts$IGIS[idx])
  idx <- is.na(gene_names) | gene_names==""
  gene_names[idx] <- paste0("unknown_gene_",which(idx))
  gene_names <- make.unique(gene_names)
  counts <- counts[,-c(1,2)]
  counts_data <- as.matrix(counts)
  storage.mode(counts_data) <- "numeric"
  if(anyNA(counts_data) || any(counts_data<0)) stop("Invalid counts in ",files[i])
  if(any(abs(counts_data-round(counts_data))>1e-8)) warning("Non-integer counts in ",files[i])
  rownames(counts_data) <- gene_names
  counts_data <- Matrix(round(counts_data),sparse=TRUE)
  seu <- CreateSeuratObject(counts=counts_data,project=sample_ids[i])
  seu$sample_id <- sample_ids[i]
  seu$sample <- sample_labels[i]
  seu$timepoint <- if(i<=2) "control" else if(i<=4) "wk12" else "wk21"
  seurat_objects[[sample_ids[i]]] <- seu
  cat(sample_labels[i],":",ncol(seu),"cells,",nrow(seu),"genes\n")
}

sce <- merge(seurat_objects[[1]],y=seurat_objects[-1],add.cell.ids=sample_ids)
sce[["percent.mt"]] <- PercentageFeatureSet(sce,pattern="^mt-")
sce[["riboRatio"]] <- PercentageFeatureSet(sce,pattern="^Rp[sl]")

p_qc <- VlnPlot(sce,features=c("nFeature_RNA","nCount_RNA","percent.mt","riboRatio"),group.by="sample",ncol=4,pt.size=0)
ggsave(file.path(output_dir,"qc_metrics_violin.pdf"),p_qc,width=12,height=4)

# 2. QC filtering
seurat_clean <- subset(sce,subset=nFeature_RNA>500 & nFeature_RNA<7000 & nCount_RNA>1000 & percent.mt<20 & riboRatio<10)

saveRDS(seurat_clean,file.path(output_dir,"QC_filtered_seurat.rds"))

# 3. Split QC-filtered cells by sample and prepare for DoubletFinder
split_seurat <- SplitObject(seurat_clean,split.by="orig.ident")
split_seurat <- split_seurat[sample_ids]

for(i in seq_along(split_seurat)){
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]],verbose=FALSE)
  split_seurat[[i]] <- FindVariableFeatures(split_seurat[[i]],nfeatures=2000,verbose=FALSE)
  split_seurat[[i]] <- ScaleData(split_seurat[[i]],vars.to.regress="percent.mt",verbose=FALSE)
  split_seurat[[i]] <- RunPCA(split_seurat[[i]],npcs=30,verbose=FALSE)
  split_seurat[[i]] <- FindNeighbors(split_seurat[[i]],reduction="pca",dims=1:30,verbose=FALSE)
  split_seurat[[i]] <- FindClusters(split_seurat[[i]],resolution=0.8,verbose=FALSE)
  split_seurat[[i]] <- RunUMAP(split_seurat[[i]],reduction="pca",dims=1:30,verbose=FALSE)
}

# 4. DoubletFinder for each sample after QC
doublet_summary <- list()
for(i in seq_along(split_seurat)){
  obj <- split_seurat[[i]]
  n_cells <- ncol(obj)
  cat("Running DoubletFinder:",sample_labels[i],"(",n_cells,"cells)\n")
  sweep_res <- paramSweep_v3(obj,PCs=1:30,sct=FALSE)
  sweep_stats <- summarizeSweep(sweep_res,GT=FALSE)
  bcmvn <- find.pK(sweep_stats)
  if(!nrow(bcmvn) || all(is.na(bcmvn$MeanBC))) stop("Could not determine pK for ",sample_labels[i])
  pK_value <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$MeanBC)]))
  homotypic_prop <- modelHomotypic(obj$seurat_clusters)
  nExp_poi <- round(n_cells*8e-6*n_cells)
  nExp_poi.adj <- round(nExp_poi*(1-homotypic_prop))
  obj <- doubletFinder_v3(obj,PCs=1:30,pN=0.25,pK=pK_value,nExp=nExp_poi.adj,reuse.pANN=FALSE,sct=FALSE)
  split_seurat[[i]] <- obj
  df_col <- tail(grep("^DF.classifications_",colnames(obj@meta.data),value=TRUE),1)
  status <- as.character(obj@meta.data[[df_col]])
  n_doublet <- sum(status=="Doublet",na.rm=TRUE)
  n_singlet <- sum(status=="Singlet",na.rm=TRUE)
  doublet_summary[[i]] <- data.frame(sample=sample_labels[i],cells_after_QC=n_cells,pK=pK_value,homotypic_prop=homotypic_prop,expected_doublets=nExp_poi,expected_doublets_adjusted=nExp_poi.adj,singlets=n_singlet,doublets=n_doublet,doublet_percent=100*n_doublet/max(n_doublet+n_singlet,1))
}

summary_df <- do.call(rbind,doublet_summary)

saveRDS(split_seurat,file.path(output_dir,"DoubletFinder_split_objects.rds"))

# 5. Plot DoubletFinder results
plots <- lapply(seq_along(split_seurat),function(i){
  obj <- split_seurat[[i]]
  df_col <- tail(grep("^DF.classifications_",colnames(obj@meta.data),value=TRUE),1)
  status <- as.character(obj@meta.data[[df_col]])
  n_doublet <- sum(status=="Doublet",na.rm=TRUE)
  n_singlet <- sum(status=="Singlet",na.rm=TRUE)
  pct <- 100*n_doublet/max(n_doublet+n_singlet,1)
  DimPlot(obj,reduction="umap",group.by=df_col,cols=c("Singlet"="skyblue","Doublet"="#E97B79"),raster=TRUE,pt.size=4,raster.dpi=c(720,720))+
    labs(title=sample_labels[i],subtitle=sprintf("Singlet: %d | Doublet: %d (%.1f%%)",n_singlet,n_doublet,pct))+
    theme(plot.title=element_text(hjust=0.5,size=8),plot.subtitle=element_text(hjust=0.5,size=8))
})

p <- wrap_plots(plots,ncol=3)






