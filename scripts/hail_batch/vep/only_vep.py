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
    # filter to biallelic loci only
    tob_wgs = tob_wgs.filter_rows(hl.len(tob_wgs.alleles) == 2)
    tob_wgs = tob_wgs.filter_rows(tob_wgs.alleles[1] != '*')
    vep = hl.vep(tob_wgs)
    vep_filename = 'gs://cpg-tob-wgs-test/kat/v0/tobwgs_v7_vep_full.mt'
    vep.write(vep_filename)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
