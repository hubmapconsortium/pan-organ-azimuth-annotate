#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from subprocess import run, check_call
from typing import Optional
import anndata
import scanpy as sc
import muon as mu

def main(
        secondary_analysis_matrix: Path,
):
    if secondary_analysis_matrix.suffix == ".h5mu":
        mudata = mu.read_h5mu(secondary_analysis_matrix)
        adata = mudata.mod["rna"]
    else:
        adata = anndata.read_h5ad(secondary_analysis_matrix)

    adata = adata[:,~adata.var.hugo_symbol.isna()]
    adata.write('secondary_analysis_hugo.h5ad')
    azimuth_annotate_command = f"annotate secondary_analysis_hugo.h5ad -fn hugo_symbol"
    check_call(azimuth_annotate_command, shell=True)
    ct_adata = anndata.read_h5ad('secondary_analysis_hugo_ANN.h5ad')
    secondary_analysis_adata = anndata.AnnData(X=adata.X, var=adata.var, obs=ct_adata.obs,
                                               obsm = ct_adata.obsm, uns=ct_adata.uns)


    sc.pl.umap(secondary_analysis_adata, color="azimuth_fine", show=False, save="umap_by_cell_type.pdf")

    if "X_spatial" in adata.obsm:
        sc.pl.scatter(secondary_analysis_adata, color="azimuth_fine", basis="spatial", show=False,
                      save="spatial_pos_by_cell_type.pdf")

    single_sample_cells = [c for c in secondary_analysis_adata.obs.azimuth_fine.unique() if
                           len(secondary_analysis_adata[secondary_analysis_adata.obs.azimuth_fine == c].obs.index) == 1]
    secondary_analysis_adata_subset = secondary_analysis_adata[~secondary_analysis_adata.obs.azimuth_fine.isin(single_sample_cells)]

    sc.tl.rank_genes_groups(secondary_analysis_adata_subset, "azimuth_fine", method="t-test",
                            key_added="rank_genes_groups_cell_type")

    sc.pl.rank_genes_groups(secondary_analysis_adata_subset, key="rank_genes_groups_cell_type", n_genes=25, sharey=False,
                            save="marker_genes_by_cell_type_t_test.pdf")

    secondary_analysis_adata.uns = secondary_analysis_adata_subset.uns


    if secondary_analysis_matrix.suffix == ".h5mu":
        mudata.mod["rna"] = secondary_analysis_adata
        mudata.write_h5mu("secondary_analysis.h5mu")
    else:
        secondary_analysis_adata.write("secondary_analysis.h5ad")

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('secondary_analysis_matrix', type=Path)
    args = p.parse_args()

    main(
        args.secondary_analysis_matrix,
    )
