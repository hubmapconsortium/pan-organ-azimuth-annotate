#!/usr/bin/env python3More actions
from argparse import ArgumentParser
from pathlib import Path
from subprocess import run, check_call
from typing import Optional
import anndata
import pandas as pd
import scanpy as sc
import squidpy as sq
import muon as mu
from plot_utils import new_plot
import matplotlib.pyplot as plt
import json

CLID_MAPPING = "/opt/pan-human-azimuth-crosswalk.csv"


def map_to_clid(adata_obs: pd.DataFrame):
    reference = pd.read_csv(CLID_MAPPING, header=10)
    obs_w_clid = adata_obs.merge(reference[['Annotation_Label', 'CL_Label', 'CL_ID']],
                                 how='left',
                                 left_on='final_level_labels',
                                 right_on='Annotation_Label')
    obs_w_clid = obs_w_clid.drop('Annotation_Label', axis='columns')
    return obs_w_clid


def main(
        secondary_analysis_matrix: Path,
        organism: str,
):
    if organism == "human":
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

        secondary_analysis_adata.obs = map_to_clid(secondary_analysis_adata.obs)
        secondary_analysis_adata.uns["pan_human_azimuth_crosswalk"] = {
            "title": "Cell type annotations for pan-human Azimuth, v1.0",
            "description": (
                "This crosswalk maps cell type annotations from pan-human Azimuth to the "
                "Cell Ontology (Version IRI: http://purl.obolibrary.org/obo/cl/releases/2025-04-10/cl.owl)."
            ),
            "url": "https://cdn.humanatlas.io/digital-objects/ctann/pan-human-azimuth/latest/assets/pan-human-azimuth-crosswalk.csv",
            "publisher": "HuBMAP",
            "creators": ["Aleix Puig-Barbe"],
            "project_lead": "Katy Börner",
            "reviewers": ["Bruce W. Herr II", "Katy Börner"],
            "processor": "HRA Digital Object Processor, v0.9.0",
            "date_published": "2025-06-15",
            "date_last_processed": "2025-06-12",
            "funders": [
                "National Institutes of Health (OT2OD033756)",
                "National Institutes of Health (OT2OD026671)"
            ],
            "license": "CC BY 4.0",
            "dashboard": "https://apps.humanatlas.io/dashboard/data"
        }

        for key in adata.obsm:
            secondary_analysis_adata.obsm[key] = adata.obsm[key]

        if 'X_umap' in secondary_analysis_adata.obsm:
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

        calculated_metadata_dict = {"annotation_tools": ["Azimuth"], "object_types": ["CL:0000000"]}
        with open('calculated_metadata.json', 'w') as f:
            json.dump(calculated_metadata_dict, f)

        cell_type_manifest_dict = {}

        for column_header in ['azimuth_broad', 'azimuth_medium', 'azimuth_fine', 'final_level_labels', 'CL_ID']:
            sub_dict = {
                val: int((secondary_analysis_adata.obs[column_header] == val).sum())
                for val in secondary_analysis_adata.obs[column_header].unique()
            }
            # Remove NaN key if it exists
            sub_dict = {k: v for k, v in sub_dict.items() if not pd.isna(k)}
            cell_type_manifest_dict[column_header] = sub_dict

        with open('cell_type_manifest.json', 'w') as f:
            json.dump(cell_type_manifest_dict, f)

    else:
        if secondary_analysis_matrix.suffix == ".h5mu":
            mudata = mu.read_h5mu(secondary_analysis_matrix)
            mudata.write_h5mu("secondary_analysis.h5mu")
        else:
            adata = anndata.read_h5ad(secondary_analysis_matrix)
            adata.write("secondary_analysis.h5ad")

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('secondary_analysis_matrix', type=Path)
    p.add_argument('organism', type=str)
    args = p.parse_args()

    main(
        args.secondary_analysis_matrix,
        args.organism,
    )
    