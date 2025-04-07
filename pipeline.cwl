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

steps:
  pan_organ_azimuth:
    run: steps/pan-organ-azimuth.cwl
    in:
      secondary_analysis_matrix: secondary_analysis_matrix
    out:
      [annotated_secondary_analysis_matrix]

