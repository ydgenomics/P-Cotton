# convert_format.R 1011
# Â© EMBL-European Bioinformatics Institute, 2023
# Yuyao Song <ysong@ebi.ac.uk>
library(optparse)
library(Seurat)
option_list <- list(
  make_option(c("-i", "--input_file"),
    type = "character", default = NULL,
    help = "Path to input file for convrting"
  ),
  make_option(c("-o", "--output_file"),
    type = "character", default = NULL,
    help = "Output file after conversion"
  ),
  make_option(c("-s", "--stype"),
    type = "character", default = NULL,
    help = "Conversion type, choose between anndata_to_seurat or seurat_to_anndata"
  )
#    make_option(c("--conda_path"),
#    type = "character", default = NULL,
#    help = "Conda for python executable to use for reticulate, important to match the prepared conda env!"
#  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if(FALSE){
opt <- list()
opt$input_file <- "/data/input/Files/ResultData/Workflow/W202409020001454/Cer_leaf__soupx_dataget/Cer_leaf__soupx.h5ad"
opt$output_file <- "/data/work/Cer_leaf__soupx.RDS"
opt$stype <- "anndata_to_seurat"
#opt$conda_path <- "~/anaconda3/envs/sceasy"
}
library(schard)
input_file <- opt$input_file
output_file <- opt$output_file
stype <- opt$stype
#conda_path <- opt$conda_path

# set sys env before loading reticulate
#Sys.setenv(RETICULATE_PYTHON=paste0(conda_path, "/bin/python3"))
#Sys.setenv(RETICULATE_PYTHON_ENV=conda_path)

#library(reticulate)
#library(sceasy)
#library(anndata)

if(stype == 'anndata_to_seurat'){

    message(paste0("from anndata to seurat, input: ", input_file))

    dt=schard::h5ad2seurat(input_file)
    saveRDS(dt,file=opt$output_file)
} else if (stype == 'seurat_to_anndata'){
    warning("using SeuratDisk and  another script")
}
