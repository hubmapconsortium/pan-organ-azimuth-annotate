#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from subprocess import run, check_call
from typing import Optional
import anndata
import scanpy

def main(
        h5ad_file: Path,
):
    azimuth_annotate_command = f"annotate {h5ad_file}"
    check_call(azimuth_annotate_command, shell=True)
    mv_command = f"mv secondary_analysis_ANN.h5ad secondary_analysis.h5ad"
    check_call(mv_command, shell=True)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('h5ad_file', type=Path)
    args = p.parse_args()

    main(
        args.h5ad_file,
    )
