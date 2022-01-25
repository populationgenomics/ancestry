#!/usr/bin/env python3


"""
Run VEP on the TOB dataset
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
    # filter to biallelic loci only
    tob_wgs = tob_wgs.filter_rows(hl.len(tob_wgs.alleles) == 2)
    tob_wgs = tob_wgs.filter_rows(tob_wgs.alleles[1] != '*')
    vep = hl.vep(tob_wgs)
    vep_path = output_path('tobwgs_vep95_GRCh38.mt')
    vep.write(vep_path)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
