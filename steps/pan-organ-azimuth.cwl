#!/usr/bin/env cwl-runner
cwlVersion: v1.1
class: CommandLineTool
requirements:
  DockerRequirement:
    dockerPull: hubmap/pan-organ-azimuth:0.1.3

inputs:
  secondary_analysis_matrix:
    type: File
    inputBinding:
      position: 0
  organism:
    type: string?
    inputBinding:
      position: 1


outputs:
  annotated_secondary_analysis_matrix:
    type: File
    outputBinding:
      glob: 'secondary_analysis.*'
  umap_plot:
    type: File?
    outputBinding:
      glob: umap_by_cell_type.pdf
  spatial_plot:
    type: File?
    outputBinding:
      glob: spatial_pos_by_cell_type.pdf
  marker_gene_plot:
    type: File?
    outputBinding:
      glob: marker_genes_by_cell_type_t_test.pdf
  neighborhood_enrichment_plot:
    type: File?
    outputBinding:
      glob: neighborhood_enrichment_by_cell_type.pdf
  calculated_metadata_file:
    type: File?
    outputBinding:
      glob: calculated_metadata.json
  cell_type_manifest:
    type: File?
    outputBinding:
      glob: cell_type_manifest.json

baseCommand: ['python3', '/opt/pan_organ_azimuth.py']
