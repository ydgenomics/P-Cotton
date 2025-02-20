# read.py concat biosample 1017
import pandas as pd
import scanpy as sc
import anndata as ad
import sys
files_txt_path = sys.argv[1]
projects_txt_path = sys.argv[2]
species=sys.argv[3]
with open(files_txt_path, 'r') as file:
    file_content = file.read().strip()
indataget = file_content.split(',')
with open(projects_txt_path, 'r') as filen:
    file_content = filen.read().strip()
projects = file_content.split(',')
print(indataget)
print(projects)
adatas={}
for i in range(len(indataget)): 
    key = projects[i]
    value = indataget[i]
    value = sc.read_h5ad(value)
    #ensure concat used by raw data
    value.X = value.layers["counts"]
    adatas[key] = value
#if len(adatas) == 1:
#    adata = list(adatas.values())[0]
#    adata.obs['biosample'] = 'biosample_value'
#elif len(adatas) > 1:
    adata = ad.concat(adatas, label="biosample")
adata.obs_names_make_unique()
print(adata.obs["biosample"].value_counts())

print(adata.obs.columns)
adata.write_h5ad(filename=species+'.h5ad',compression="gzip")
