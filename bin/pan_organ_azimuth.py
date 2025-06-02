#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from subprocess import run, check_call
from typing import Optional
import anndata
import scanpy as sc
import squidpy as sq
import muon as mu
from plot_utils import new_plot
import matplotlib.pyplot as plt

def main(
        secondary_analysis_matrix: Path,
):
    if secondary_analysis_matrix.suffix == ".h5mu":
        mudata = mu.read_h5mu(secondary_analysis_matrix)
        adata = mudata.mod["rna"]
    else:
        adata = anndata.read_h5ad(secondary_analysis_matrix)

    if "unscaled" in adata.layers:
        adata.X = adata.layers["unscaled"]

    adata = adata[:,~adata.var.hugo_symbol.isna()]

    adata.write('secondary_analysis_hugo.h5ad')
    azimuth_annotate_command = f"annotate secondary_analysis_hugo.h5ad -fn hugo_symbol"
    check_call(azimuth_annotate_command, shell=True)
    ct_adata = anndata.read_h5ad('secondary_analysis_hugo_ANN.h5ad')
    secondary_analysis_adata = anndata.AnnData(X=adata.X, var=adata.var, obs=ct_adata.obs,
                                               obsm = ct_adata.obsm, uns=ct_adata.uns)
    for key in adata.obsm:
        secondary_analysis_adata.obsm[key] = adata.obsm[key]

    with new_plot():
        sc.pl.umap(secondary_analysis_adata, color="final_level_labels", show=False)
        plt.savefig("umap_by_cell_type.pdf", bbox_inches="tight")

    if "X_spatial" in adata.obsm:
        if "spatial" not in adata.obsm:
            secondary_analysis_adata.obsm["spatial"] = secondary_analysis_adata.obsm["X_spatial"]

        with new_plot():
            sc.pl.scatter(secondary_analysis_adata, color="final_level_labels", basis="spatial", show=False)
            plt.savefig("spatial_pos_by_cell_type.pdf", bbox_inches="tight")

        sq.gr.spatial_neighbors(secondary_analysis_adata)
        secondary_analysis_adata_subset = secondary_analysis_adata[~secondary_analysis_adata.obs.final_level_labels.isna()]

        sq.gr.nhood_enrichment(secondary_analysis_adata_subset, cluster_key="final_level_labels")

        with new_plot():
            sq.pl.nhood_enrichment(secondary_analysis_adata_subset, cluster_key="final_level_labels")
            plt.savefig("neighborhood_enrichment_by_cell_type.pdf", bbox_inches="tight")

#    single_sample_cells = [c for c in secondary_analysis_adata.obs.azimuth_fine.unique() if
#                           len(secondary_analysis_adata[secondary_analysis_adata.obs.azimuth_fine == c].obs.index) == 1]
#    secondary_analysis_adata_subset = secondary_analysis_adata[~secondary_analysis_adata.obs.azimuth_fine.isin(single_sample_cells)]

#    sc.tl.rank_genes_groups(secondary_analysis_adata_subset, "azimuth_fine", method="t-test",
#                            key_added="rank_genes_groups_cell_type")

#    with new_plot():
#        sc.pl.rank_genes_groups(secondary_analysis_adata_subset, key="rank_genes_groups_cell_type",
##                                n_genes=25, sharey=False)
 #       plt.savefig("marker_genes_by_cell_type_t_test.pdf", bbox_inches="tight")

#    secondary_analysis_adata.uns = secondary_analysis_adata_subset.uns


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

