# SCTransform.4methods_integration.R 250306
# SCTransform.harmony_integration.R 250225
# https://satijalab.org/seurat/articles/seurat5_integration
# Interesting thing is written for V5.20 'split()' and 'IntegrateLayers'
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
  make_option(c("-r", "--resolution_set"),
    type = "double", default = NULL,
    help = "Set the resolution for clustering"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input_rds)){
  opt$input_rds <- "/data/work/convert/Cer_test_convert.rds"
  opt$out_rds <- "/data/work/sct_harmony/Cer_test_convert_SCTransform.harmony_integrated.rds" # nolint
  opt$out_UMAP <- "/data/work/sct_harmony/Cer_test_convert_SCTransform.harmony_integrated_UMAP.pdf"  # nolint
  opt$batch_key <- "biosample"
  opt$sample_key <- "sample"
  opt$ resolution_set <- 1.0
}

input_rds <- opt$input_rds
out_rds <- opt$out_rds
out_UMAP <- opt$out_UMAP
batch_key <- opt$batch_key
sample_key <- opt$sample_key
resolution_set <- opt$ resolution_set

obj <- readRDS(input_rds)
#obj <- subset(obj, nFeature_RNA > 1000)

obj[["RNA"]] <- split(obj[["RNA"]], f = obj$biosample)

# run sctransform
obj <- SCTransform(obj, vst.flavor = "v2")
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)

# Perform streamlined (one-line) integrative analysis
obj <- IntegrateLayers(
  object = obj, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "cca",
  assay = "SCT", verbose = FALSE
)
obj <- IntegrateLayers(
  object = obj, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "rpca",
  assay = "SCT", verbose = FALSE
)
obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  assay = "SCT", verbose = FALSE
)
obj <- IntegrateLayers(
  object = obj, method = FastMNNIntegration,
  new.reduction = "mnn",
  assay = "SCT", verbose = FALSE
)

out_vln <- gsub("_UMAP.pdf", "_VlnPlot.pdf", out_UMAP)
pdf(out_vln)
VlnPlot(obj, features = c("nCount_RNA", "nFeature_RNA"), group.by= batch_key)
dev.off()

out_UMAP <- gsub(".pdf", "", out_UMAP) # nolint

integrations <- c("cca", "rpca", "harmony", "mnn") # nolint
cluster_names <- c("cca_clusters", "rpca_clusters", "harmony_clusters", "mnn_clusters") # nolint

for (i in seq_along(integrations)) {
  obj <- FindNeighbors(obj, reduction = integrations[i], dims = 1:30)
  obj <- FindClusters(obj, resolution = resolution_set, cluster.name = cluster_names[i]) # nolint
  obj <- RunUMAP(obj, reduction = integrations[i], dims = 1:30, reduction.name = paste0("umap_", integrations[i])) # nolint
  pdf(paste0(out_UMAP, "_", integrations[i], ".pdf"))
  DimPlot(obj, reduction = paste0("umap_", integrations[i]), split.by = batch_key) # nolint
  DimPlot(obj, reduction = paste0("umap_", integrations[i]), group.by = batch_key, shuffle = TRUE, label = TRUE) # nolint
  DimPlot(obj, reduction = paste0("umap_", integrations[i]), group.by = sample_key, shuffle = TRUE, label = TRUE) # nolint
  DimPlot(obj, reduction = paste0("umap_", integrations[i]), group.by = cluster_names[i], shuffle = TRUE, label = TRUE) # nolint
  dev.off()
}

DefaultAssay(obj) <- "RNA"
obj [["RNA"]] <- JoinLayers(obj [["RNA"]])
obj[["RNA"]] <- as(obj[["RNA"]], "Assay") # Assay RNA changing from Assay5 to Assay # nolint

saveRDS(obj, file = out_rds)
str(obj)
