#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.1
requirements:
  ScatterFeatureRequirement: {}

inputs:
  secondary_analysis_matrix:
    type: File

outputs:
  annotated_secondary_analysis_matrix:
    type: File
    outputSource: pan_organ_azimuth/annotated_secondary_analysis_matrix
  umap_plot:
    outputSource: pan_organ_azimuth/umap_plot
    type: File
    label: "UMAP dimensionality reduction plot colored by cell type"
  spatial_plot:
    outputSource: pan_organ_azimuth/spatial_plot
    type: File?
    label: "Spatial plot colored by cell type"
  marker_gene_plot_t_test:
    outputSource: pan_organ_azimuth/marker_gene_plot
    type: File
    label: "Cell type marker genes, t-test"
  neighborhood_enrichment_plot:
    outputSource: pan_organ_azimuth/neighborhood_enrichment_plot
    type: File?

steps:
  pan_organ_azimuth:
    run: steps/pan-organ-azimuth.cwl
    in:
      secondary_analysis_matrix: secondary_analysis_matrix
    out:
      [annotated_secondary_analysis_matrix, umap_plot, spatial_plot, marker_gene_plot, neighborhood_enrichment_plot]

