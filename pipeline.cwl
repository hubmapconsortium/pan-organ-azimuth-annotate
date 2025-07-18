#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.1
requirements:
  ScatterFeatureRequirement: {}

inputs:
  secondary_analysis_matrix:
    type: File
  organism:
    type: string?
    default: "human"

outputs:
  annotated_secondary_analysis_matrix:
    type: File
    outputSource: pan_organ_azimuth/annotated_secondary_analysis_matrix
  umap_plot:
    outputSource: pan_organ_azimuth/umap_plot
    type: File?
    label: "UMAP dimensionality reduction plot colored by cell type"
  spatial_plot:
    outputSource: pan_organ_azimuth/spatial_plot
    type: File?
    label: "Spatial plot colored by cell type"
  marker_gene_plot_t_test:
    outputSource: pan_organ_azimuth/marker_gene_plot
    type: File?
    label: "Cell type marker genes, t-test"
  neighborhood_enrichment_plot:
    outputSource: pan_organ_azimuth/neighborhood_enrichment_plot
    type: File?
  calculated_metadata_file:
    outputSource: pan_organ_azimuth/calculated_metadata_file
    type: File?
  cell_type_manifest:
    outputSource: pan_organ_azimuth/cell_type_manifest
    type: File?

steps:
  pan_organ_azimuth:
    run: steps/pan-organ-azimuth.cwl
    in:
      secondary_analysis_matrix: secondary_analysis_matrix
      organism: organism
    out:
      [annotated_secondary_analysis_matrix, umap_plot, spatial_plot, marker_gene_plot, neighborhood_enrichment_plot,
      calculated_metadata_file, cell_type_manifest]

