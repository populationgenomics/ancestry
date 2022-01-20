"""Calculate ld using the ld_matrix function"""

import hail as hl
import pandas as pd
from hail.linalg import BlockMatrix

BLOCK_MATRIX = 'gs://cpg-tob-wgs-test/kat/v0/ld_matrix_1M_100k.bm/'


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    bm = BlockMatrix.read('gs://cpg-tob-wgs-test/kat/v0/ld_matrix_1M_100k.bm/')
    bm = bm.to_numpy()
    bm = pd.DataFrame(bm)
    # get number of non-zero, non na values across rows
    nonzero_values = bm.fillna(0).astype(bool).sum(axis=1)
    nonzero_values.to_csv(
        'gs://cpg-tob-wgs-test/kat/v1/nonzero_nona_values_1M_radius_100k.csv'
    )
    # get number of positive values (including 1's) across rows
    positive_values = bm.fillna(0).gt(0).sum(axis=1)
    positive_values.to_csv(
        'gs://cpg-tob-wgs-test/kat/v1/positive_values_1M_radius_100k.csv'
    )
    # get number of negative values across rows
    negative_values = bm.fillna(0).lt(0).sum(axis=1)
    negative_values.to_csv(
        'gs://cpg-tob-wgs-test/kat/v1/negative_values_1M_radius_100k.csv'
    )


if __name__ == '__main__':
    query()
