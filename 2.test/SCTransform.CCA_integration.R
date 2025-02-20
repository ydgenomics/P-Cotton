# SCTransform.CCA_integration.R 1010
# https://satijalab.org/seurat/archive/v4.3/sctransform_v2_vignette
library(Seurat) # make sure you are running SeuratV5
options(Seurat.object.assay.version = 'v5')
library(SeuratData)
library(patchwork)
library(optparse)
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

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input_rds)){
  opt$input_rds <- "/data/work/convert/Cer_test_convert.rds"
}
if (is.null(opt$out_rds)){
  opt$out_rds <- "/data/work/sct/Cer_test_convert_SCTransform.CCA_integrated.rds"
}
if (is.null(opt$out_UMAP)){
  opt$out_UMAP <- "/data/work/sct/Cer_test_convert_SCTransform.CCA_integrated_UMAP.pdf"
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

#pre-processs
obj <- readRDS(input_rds)
obj.list <- SplitObject(obj, split.by = batch_key)
obj.list.transformed <- list()
batch_keys <- unique(obj@meta.data$biosample)
for (biosample in batch_keys) {
  current.obj <- obj.list[[biosample]]
  transformed.obj <- SCTransform(current.obj, vst.flavor = "v2", verbose = FALSE) %>% RunPCA(npcs = 30, verbose = FALSE)
  obj.list.transformed[[biosample]] <- transformed.obj
}
obj.list <- obj.list.transformed
obj.list


#Perform integration 
#features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
features <- SelectIntegrationFeatures(object.list = obj.list)
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features)
#method=cca
obj.anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT", anchor.features = features, reduction = "cca")
obj <- IntegrateData(anchorset = obj.anchors, normalization.method = "SCT")


#Perform an integrated analysis
obj <- RunPCA(obj, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- FindClusters(obj, resolution = resolution_set, cluster.name = "celltype")

#save rds
saveRDS(obj, file = out_rds)
#
str(obj)
#visual
pdf(out_UMAP)
DimPlot(obj, reduction = "umap", split.by = batch_key)
DimPlot(obj, reduction = "umap", group.by = batch_key, label = TRUE, shuffle = TRUE)
DimPlot(obj, reduction = "umap", group.by = sample_key, label = TRUE, shuffle = TRUE)
DimPlot(obj, reduction = "umap", group.by = cluster_key, label = TRUE, shuffle = TRUE)
VlnPlot(obj, features = c("nCount_SCT", "nFeature_SCT"), group.by= batch_key)
dev.off()

