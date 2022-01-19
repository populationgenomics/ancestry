#!/usr/bin/env python3


"""
python script to try and run vep
"""


import click
import hail as hl

TOB_WGS = 'gs://cpg-tob-wgs-test/mt/v7.mt/'


@click.command()
def main():
    """
    just runs vep
    """

    hl.init(default_reference='GRCh38')

    tob_wgs = hl.read_matrix_table(TOB_WGS)
    hl.vep(tob_wgs)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
