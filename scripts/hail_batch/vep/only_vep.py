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
    Run vep using run_vep.py wrapper
    """

    hl.init(default_reference='GRCh38')

    tob_wgs = hl.read_matrix_table(TOB_WGS)
    tob_wgs = tob_wgs.head(10)
    vep = hl.vep(tob_wgs)
    vep_filename = 'gs://cpg-tob-wgs-test/kat/v0/tobwgs_v7_vep.mt'
    vep.write(vep_filename)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
