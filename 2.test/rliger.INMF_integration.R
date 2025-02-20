# rliger.INMF_integration.R 1011
# https://welch-lab.github.io/liger/articles/Integrating_multi_scRNA_data.html#r-session-info
# https://github.com/Papatheodorou-Group/BENGAL/blob/main/bin/rliger_integration_UINMF_multiple_species.R
library(optparse)
library(rliger)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
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
  opt$out_rds <- "/data/work/liger/Cer_test_convert_rliger.INMF_integrated.rds"
}
if (is.null(opt$out_UMAP)){
  opt$out_UMAP <- "/data/work/liger/Cer_test_convert_rliger.INMF_integrated_UMAP.pdf"
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
#
obj <- readRDS(input_rds)
obj <- obj %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData(split.by = batch_key, do.center = FALSE)
# LIGER
obj <- RunOptimizeALS(obj, k = 30, lambda = 5, split.by = batch_key)
obj <- RunQuantileNorm(obj, split.by = batch_key)
names(obj@reductions)
#obj <- FindNeighbors(obj, reduction = "iNMF_raw", k.param = 10, dims = 1:30)
obj <- FindNeighbors(obj, reduction = "iNMF", k.param = 10, dims = 1:30)
obj <- FindClusters(obj, resolution = resolution_set, cluster.name = "celltype")
# Dimensional reduction and plotting
#obj <- RunUMAP(obj, dims = 1:ncol(obj[["iNMF_raw"]]), reduction = "iNMF_raw", n_neighbors = 15L,  min_dist = 0.3)
obj <- RunUMAP(obj, dims = 1:ncol(obj[["iNMF"]]), reduction = "iNMF", n_neighbors = 15L)
#obj <- RunUMAP(obj, reduction = "iNMF", n_neighbors = 30, min_dist = 0.3)
#for below scib_test
obj@reductions$pca <- obj@reductions$iNMF
names(obj@reductions)
# have to convert all factor to character, or when later converting to h5ad, the factors will be numbers
i <- sapply(obj@meta.data, is.factor)
obj@meta.data[i] <- lapply(obj@meta.data[i], as.character)
#
obj
# iNMF embedding will be in .obsm['iNMF']
saveRDS(obj, file = out_rds)
pdf(out_UMAP)
DimPlot(obj, reduction = "umap", split.by = batch_key)
DimPlot(obj, reduction = "umap", group.by = batch_key, shuffle = TRUE, label = TRUE)
DimPlot(obj, reduction = "umap", group.by = sample_key, shuffle = TRUE, label = TRUE)
DimPlot(obj, reduction = "umap", group.by = cluster_key, shuffle = TRUE, label = TRUE)
VlnPlot(obj, features = c("nCount_RNA", "nFeature_RNA"), group.by= batch_key)
dev.off()

