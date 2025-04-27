# BBKNNR_integration.R 1010
# /opt/conda/bin/R
#https://mp.weixin.qq.com/s/ZkY8R3yZEEsIuV8lDIAdlA
library(bbknnR)
library(Seurat)
library(ggplot2)
library(dplyr)
library(SeuratData)
library(patchwork)
library(optparse)
library(magrittr)

option_list <- list(
  make_option(c("-i", "--input_rds"),
    type = "character", default = NULL,
    help = "Path to input preprocessed rds file"
  ),
  make_option(c("-o", "--out_rds"),
    type = "character", default = NULL,
    help = "integrated rds file"
  ),
  make_option(c("-p", "--out_UMAP"),
    type = "character", default = NULL,
    help = "Output UMAP after integration"
  ),
  make_option(c("-b", "--batch_key"),
    type = "character", default = NULL,
    help = "Batch key identifier to integrate"
  ),
  make_option(c("-s", "--sample_key"),
    type = "character", default = NULL,
    help = "Sample key identifier"
  ),
  make_option(c("-c", "--cluster_key"),
    type = "character", default = NULL,
    help = "Cluster key for UMAP plotting"
  ),
  make_option(c("-r", "--resolution_set"),
    type = "double", default = NULL,
    help = "Set the resolution for clustering"
  )
)

# parse input
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input_rds)){
  opt$input_rds <- "/data/work/convert/Cer_test_convert.rds"
}
if (is.null(opt$out_rds)){
  opt$out_rds <- "/data/work/bbknn/Cer_BBKNNR_integrated.rds"
}
if (is.null(opt$out_UMAP)){
  opt$out_UMAP <- "/data/work/bbknn/Cer_BBKNNR_integrated_UMAP.pdf"
}
if (is.null(opt$batch_key)){
  opt$batch_key <- "biosample"
}
if (is.null(opt$sample_key)){
  opt$sample_key <- "sample"
}
if (is.null(opt$cluster_key)){
  opt$cluster_key <- "celltype"
}
if (is.null(opt$resolution_set)){
  opt$ resolution_set <- 1.0
}
#
input_rds <- opt$input_rds
out_rds <- opt$out_rds
out_UMAP <- opt$out_UMAP
batch_key <- opt$batch_key
sample_key <- opt$sample_key
cluster_key <- opt$cluster_key
resolution_set <- opt$ resolution_set

#### 1.load dataset
obj <- readRDS(input_rds)
#obj[["RNA"]] <- split(obj[["RNA"]], f = obj$biosample)
#### 2.normalize/HVG/scale/pca
#obj <- NormalizeData(obj) 
#obj <- FindVariableFeatures(obj, selection.method = "vst")
#obj <- ScaleData(obj)
#obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
obj <- NormalizeData(obj) %>% FindVariableFeatures(selection.method = "vst") %>% ScaleData() %>% RunPCA(npcs = 50, verbose = FALSE)

#### 3. bbknn
obj <- RunBBKNN(obj, reduction = "pca", run_TSNE = FALSE, batch_key = batch_key)
#obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
#obj <- FindNeighbors(obj, reduction = "pca", k.param = 10, dims = 1:30) 
#obj <- FindClusters(obj, resolution = resolution_set, cluster.name = "integrated_cluster", algorithm = 1, graph.name="bbknn")
#obj <- FindClusters(obj, resolution = resolution_set, cluster.name = "integrated_cluster")

obj <- FindNeighbors(obj, reduction = "pca", k.param = 10, dims = 1:30) %>%
  FindClusters(resolution = resolution_set, algorithm = 1, graph.name = "bbknn", cluster.name = "celltype") %>%
  identity()
unique(obj$celltype)
obj
#
#obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, reduction.name = "umap")
saveRDS(obj, file = out_rds)

#visual
pdf(out_UMAP)
DimPlot(obj, reduction = "umap", split.by = batch_key)
DimPlot(obj, reduction = "umap", group.by = batch_key, shuffle = TRUE, label = TRUE)
DimPlot(obj, reduction = "umap", group.by = sample_key, shuffle = TRUE, label = TRUE)
DimPlot(obj, reduction = "umap", group.by = cluster_key, shuffle = TRUE, label = TRUE)
VlnPlot(obj, features = c("nCount_RNA", "nFeature_RNA"), group.by= batch_key)
dev.off()
#VlnPlot(obj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), group.by= "biosample")