#!/usr/bin/env python3


"""
Run VEP on the hail mt
"""


import click
import hail as hl
from analysis_runner import output_path


@click.command()
@click.option('--input', required=True, help='Hail matrix table to run VEP on')
def main(input: str):
    """
    Run vep using main.py wrapper
    """

    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(input)
    # filter to biallelic loci only
    mt = mt.filter_rows(hl.len(mt.alleles) == 2)
    mt = mt.filter_rows(mt.alleles[1] != '*')
    vep = hl.vep(mt)
    vep_path = output_path('vep95_GRCh38.mt')
    vep.write(vep_path)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
