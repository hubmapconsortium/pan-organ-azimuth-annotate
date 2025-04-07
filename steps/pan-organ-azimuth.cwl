#!/usr/bin/env cwl-runner
cwlVersion: v1.1
class: CommandLineTool
requirements:
  DockerRequirement:
    dockerPull: hubmap/pan-organ-azimuth:latest

inputs:
  secondary_analysis_matrix:
    type: File
    inputBinding:
      position: 0

outputs:
  annotated_secondary_analysis_matrix:
    type: File
    outputBinding:
      glob: 'secondary_analysis.*'

baseCommand: ['python3', '/opt/pan_organ_azimuth.py']
