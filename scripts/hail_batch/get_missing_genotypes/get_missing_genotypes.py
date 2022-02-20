"""Get number of variants with a MAF > 0.05"""

import hail as hl
from analysis_runner import bucket_path

TOB_WGS = bucket_path('mt/v5.1.mt/')


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    tob_wgs = hl.read_matrix_table(TOB_WGS)
    # densify the table and split multiallelics
    tob_wgs = hl.experimental.densify(tob_wgs)
    tob_wgs = hl.split_multi_hts(tob_wgs)
    # only keep biallelic and multiallelic loci
    tob_wgs = tob_wgs.filter_rows(hl.len(tob_wgs.alleles) >= 2)
    # save file
    tob_wgs.write('gs://cpg-tob-wgs-test/kat/v0/tob_wgs_densified_filtered.mt')


if __name__ == '__main__':
    query()
