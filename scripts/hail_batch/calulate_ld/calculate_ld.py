"""Calculate ld using the ld_matrix function"""

import hail as hl
from hail.linalg import BlockMatrix

TOB_WGS = 'gs://cpg-tob-wgs-test/mt/v7.mt/'


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    tob_wgs = hl.read_matrix_table(TOB_WGS)
    tob_wgs = hl.experimental.densify(tob_wgs)
    tob_wgs = tob_wgs.head(10000)
    ld = hl.ld_matrix(tob_wgs.GT.n_alt_alleles(), tob_wgs.locus, radius=2e6)
    BlockMatrix.write(ld, 'gs://cpg-tob-wgs-test/kat/v0/ld_matrix.bm')


if __name__ == '__main__':
    query()
