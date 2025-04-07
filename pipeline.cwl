#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.1
requirements:
  ScatterFeatureRequirement: {}

inputs:
  h5ad_file:
    type: File

outputs:
  annotated_h5ad_file:
    type: File
    outputSource: pan_organ_azimuth/annotated_h5ad_file

steps:
  pan_organ_azimuth:
    run: steps/pan-organ-azimuth.cwl
    in:
      h5ad_file: h5ad_file
    out:
      [annotated_h5ad_file]

