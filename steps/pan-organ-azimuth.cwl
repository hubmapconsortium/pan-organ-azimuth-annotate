#!/usr/bin/env cwl-runner
cwlVersion: v1.1
class: CommandLineTool
requirements:
  DockerRequirement:
    dockerPull: hubmap/pan-organ-azimuth:latest

inputs:
  h5ad_file:
    type: File
    inputBinding:
      position: 0

outputs:
  annotated_h5ad_file:
    type: File
    outputBinding:
      glob: 'secondary_analysis.h5ad'

baseCommand: ['python3', '/opt/pan-organ-azimuth.py']
