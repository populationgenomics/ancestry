"""Calculate ld using the ld_matrix function"""

import hail as hl
from analysis_runner import bucket_path

TOB_WGS = bucket_path('mt/v7.mt/')


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    tob_wgs = hl.read_matrix_table(TOB_WGS)
    tob_wgs = hl.experimental.densify(tob_wgs)
    # filter out constant variants
    tob_wgs = tob_wgs.filter_rows(hl.len(tob_wgs.alleles) == 2)
    tob_wgs = tob_wgs.head(50000)
    ld = hl.ld_matrix(tob_wgs.GT.n_alt_alleles(), tob_wgs.locus, radius=2e6)
    # save block matrix
    ld.write('gs://cpg-tob-wgs-test/kat/v0/ld_matrix_2M_50k.bm')


if __name__ == '__main__':
    query()


if __name__ == '__main__':
    query()
