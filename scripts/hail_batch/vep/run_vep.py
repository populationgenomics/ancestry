#!/usr/bin/env python3


"""
Run VEP on the hail mt
"""


import click
import hail as hl
from cpg_utils.hail_batch import output_path


@click.command()
@click.option('--mt', required=True, help='Hail matrix table to run VEP on')
@click.option('--vep-version', help='Version of VEP', default='104.3')
def main(mt: str, vep_version: str):
    """
    Run vep using main.py wrapper
    """

    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(mt)
    # filter to biallelic loci only
    mt = mt.filter_rows(hl.len(mt.alleles) == 2)
    mt = mt.filter_rows(mt.alleles[1] != '*')
    vep = hl.vep(mt, config='file:///vep_data/vep-gcloud.json')
    vep_path = output_path(f'vep{vep_version}_GRCh38.mt')
    vep.write(vep_path)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
