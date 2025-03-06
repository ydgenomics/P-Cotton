# _integration.R 250223
# https://satijalab.org/seurat/articles/seurat5_integration
# Interesting thing is written for V5.20 'split()' and 'IntegrateLayers'

library(Seurat)
library(SeuratData)
library(SeuratWrappers)
# library(Azimuth)
library(ggplot2)
library(optparse)
library(patchwork)

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
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input_rds)){
  opt$input_rds <- "/data/work/convert/Cer_test_convert.rds"
}
if (is.null(opt$out_rds)){
  opt$out_rds <- "/data/work/sct_harmony/Cer_test_convert_SCTransform.harmony_integrated.rds"
}
if (is.null(opt$out_UMAP)){
  opt$out_UMAP <- "/data/work/sct_harmony/Cer_test_convert_SCTransform.harmony_integrated_UMAP.pdf"
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
if (is.null(opt$ resolution_set)){
  opt$ resolution_set <- 1.0
}

input_rds <- opt$input_rds
out_rds <- opt$out_rds
out_UMAP <- opt$out_UMAP
batch_key <- opt$batch_key
sample_key <- opt$sample_key
cluster_key <- opt$cluster_key
resolution_set <- opt$ resolution_set

obj <- readRDS(input_rds)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$Method)
obj

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)

#obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
#obj <- FindClusters(obj, resolution = resolution_set, cluster.name = "unintegrated_clusters")

#obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth
#DimPlot(obj, reduction = "umap.unintegrated", group.by = c("Method", "predicted.celltype.l2"))

#obj <- IntegrateLayers(
#  object = obj, method = CCAIntegration,
#  orig.reduction = "pca", new.reduction = 'integrated.cca',
#  verbose = FALSE)

obj <- IntegrateLayers(
  object = obj, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = 'integrated.rpca',
  verbose = FALSE)

#obj <- IntegrateLayers(
#  object = obj, method = HarmonyIntegration,
#  orig.reduction = "pca", new.reduction = 'harmony',
#  verbose = FALSE)

obj <- IntegrateLayers(
  object = obj, method = FastMNNIntegration,
  new.reduction = 'integrated.mnn',
  verbose = FALSE)

obj <- FindNeighbors(obj, reduction = 'integrated.rpca', dims = 1:30)
obj <- FindClusters(obj,resolution = 2, cluster.name = 'rpca_clusters')
obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:30, reduction.name = 'umap.cca')
p1 <- DimPlot(
  obj, reduction = "umap.cca",
  group.by = c("Method", "predicted.celltype.l2", "cca_clusters"),
  combine = FALSE, label.size = 2) 

obj <- FindNeighbors(obj, reduction = 'integrated.scvi', dims = 1:30)
obj <- FindClusters(obj,resolution = 2, cluster.name = 'scvi_clusters')
obj <- RunUMAP(obj, reduction = "integrated.scvi", dims = 1:30, reduction.name = 'umap.scvi')
p2 <- DimPlot(
  obj, reduction = "umap.scvi",
  group.by = c("Method", "predicted.celltype.l2", "scvi_clusters"),
  combine = FALSE, label.size = 2)

wrap_plots(c(p1, p2), ncol = 2, byrow = F)

p1 <- VlnPlot(
  obj, features = "rna_CD8A", group.by = 'unintegrated_clusters'
) + NoLegend() + ggtitle("CD8A - Unintegrated Clusters")
p2 <- VlnPlot(
  obj, "rna_CD8A", group.by = 'cca_clusters'
) + NoLegend() + ggtitle("CD8A - CCA Clusters")
p3 <- VlnPlot(
  obj, "rna_CD8A", group.by = 'scvi_clusters'
) + NoLegend() + ggtitle("CD8A - scVI Clusters")
p1 | p2 | p3

obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:30, reduction.name = 'umap.rpca')
p4 <- DimPlot(obj, reduction="umap.unintegrated", group.by=c("cca_clusters"))
p5 <- DimPlot(obj, reduction="umap.rpca", group.by=c("cca_clusters"))
p6 <- DimPlot(obj, reduction="umap.scvi", group.by=c("cca_clusters"))
p4 | p5 | p6

obj <- JoinLayers(obj)
obj

saveRDS(obj, file=out_rds)

