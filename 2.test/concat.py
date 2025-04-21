# Title: concat.py
# Date: 2025-04-21
# Fix: add `join` key
import pandas as pd
import scanpy as sc
import anndata as ad
import sys
files_txt_path = sys.argv[1]
projects_txt_path = sys.argv[2]
species=sys.argv[3]
group_key=sys.argv[4] #"biosample"
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
    value.X = value.layers["counts"] #ensure concat used by raw data
    adatas[key] = value
    adata = ad.concat(adatas, label=group_key, join="outer")
adata.obs_names_make_unique()
print(adata.obs[group_key].value_counts())
print(adata.obs.columns)
adata.write_h5ad(filename=species+'.h5ad',compression="gzip")
