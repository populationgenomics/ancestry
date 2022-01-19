#!/usr/bin/env python3


"""
python script to try and run vep
"""


import click
import hail as hl
from analysis_runner import bucket_path, output_path

TOB_WGS = bucket_path('mt/v7.mt/')


@click.command()
def main():
    """
    Run vep using run_vep.py wrapper
    """

    hl.init(default_reference='GRCh38')

    tob_wgs = hl.read_matrix_table(TOB_WGS)
    vep = hl.vep(tob_wgs)
    vep_filename = output_path(f'tobwgs_v7_vep.mt')
    vep.write(vep_filename)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
