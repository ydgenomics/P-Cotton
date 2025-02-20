# SCTransform.harmony_integration.R 1106
library(Seurat) # make sure you are running SeuratV5
options(Seurat.object.assay.version = 'v5')
library(SeuratData)
library(patchwork)
library(optparse)
library(ggplot2)

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
#obj <- subset(obj, nFeature_RNA > 1000)

obj[["RNA"]] <- split(obj[["RNA"]], f = obj$biosample)

# run sctransform
obj <- SCTransform(obj, vst.flavor = "v2")
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)

# one-liner to run Integration
obj <- IntegrateLayers(object = obj, method = HarmonyIntegration,
                       orig.reduction = "pca", new.reduction = 'harmony',
                       assay = "SCT", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
# obj <- FindClusters(obj, resolution = 2, cluster.name = "harmony_clusters")
obj <- FindClusters(obj, resolution = resolution_set, cluster.name = "celltype")

#colnames(obj@meta.data)[colnames(obj@meta.data) == "_index"] <- "X_index"
#
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap")

DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)

saveRDS(obj, file = out_rds)
#
unique(obj$celltype)
str(obj)
#
pdf(out_UMAP)
DimPlot(obj, reduction = "umap", split.by = batch_key)
DimPlot(obj, reduction = "umap", group.by = batch_key, shuffle = TRUE, label = TRUE)
DimPlot(obj, reduction = "umap", group.by = sample_key, shuffle = TRUE, label = TRUE)
DimPlot(obj, reduction = "umap", group.by = cluster_key, shuffle = TRUE, label = TRUE)
VlnPlot(obj, features = c("nCount_RNA", "nFeature_RNA"), group.by= batch_key)
dev.off()
