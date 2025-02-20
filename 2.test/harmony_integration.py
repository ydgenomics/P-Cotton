## harmony_integration.py 1025
# https://github.com/Papatheodorou-Group/BENGAL/blob/main/bin/harmony_integration.py
import click
import matplotlib.pyplot as plt
import scanpy as sc
import pandas
import harmonypy
import leidenalg
from matplotlib.backends.backend_pdf import PdfPages

@click.command()
@click.argument("input_h5ad", type=click.Path(exists=True))
@click.argument("out_h5ad", type=click.Path(exists=False), default=None)
@click.argument("out_umap", type=click.Path(exists=False), default=None)
@click.option('--batch_key', type=str, default=None, help="Batch key in identifying HVG and harmony integration")
@click.option('--sample_key', type=str, default=None, help="Sample key to distinguish sample")
@click.option('--cluster_key', type=str, default=None, help="Cluster key in samples")
@click.option('--resolution_set', type=float, default=None, help="set for resolution, is float")

def run_harmony(input_h5ad, out_h5ad, out_umap, batch_key, sample_key, cluster_key, resolution_set):
    click.echo('Start harmony integration')
    sc.set_figure_params(dpi_save=300, frameon=False, figsize=(10, 6))
    #input_h5ad="/data/work/read/Cer_test.h5ad"
    adata = sc.read_h5ad(input_h5ad)
    adata.var_names_make_unique()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    click.echo("HVG")
    sc.pp.highly_variable_genes(adata, batch_key=batch_key)
    #sc.pp.scale(adata, max_value=10)
    sc.pp.scale(adata)
    sc.tl.pca(adata)
    # sc.tl.pca(input_ad, svd_solver="arpack", mask_var="highly_variable")
    sc.pp.neighbors(adata, use_rep='X_pca', n_neighbors=15, n_pcs=40)
    sc.tl.umap(adata, min_dist=0.3) ## to match min_dist in seurat
    #adata.obsm['X_umapraw'] = adata.obsm['X_umap']
    click.echo("Harmony")
    sc.external.pp.harmony_integrate(adata, key=batch_key, basis = 'X_pca')
    #sc.pp.neighbors(adata, use_rep='X_pca_harmony', key_added = 'harmony', n_neighbors=15, n_pcs=40)
    sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_neighbors=15, n_pcs=40)
    sc.tl.leiden(adata, resolution=resolution_set, key_added='celltype', flavor='igraph', n_iterations=2, directed=False) 
    adata.obs['celltype'].unique()
    sc.tl.umap(adata, neighbors_key = 'neighbors') ## to match min_dist in seurat
    with PdfPages(out_umap) as pdf:
        sc.pl.umap(adata, color=[batch_key, sample_key, cluster_key], ncols=1)
        plt.savefig(pdf, format='pdf', dpi=300, bbox_inches='tight')
        plt.close()
        sc.pl.violin(adata, keys=['total_counts'], log=False, groupby='biosample')
        plt.savefig(pdf, format='pdf', dpi=300, bbox_inches='tight')
        plt.close()
        sc.pl.violin(adata, keys=['n_genes'], log=False, groupby='biosample')
        plt.savefig(pdf, format='pdf', dpi=300, bbox_inches='tight')
        plt.close()
    #adata.obsm['X_umapharmony'] = adata.obsm['X_umap']
    #click.echo("scvi integrated adata structure")
    #
    adata
    click.echo("Save output")
    adata.write(filename=out_h5ad,compression="gzip")
    click.echo("Done harmony")


if __name__ == '__main__':
    run_harmony()
