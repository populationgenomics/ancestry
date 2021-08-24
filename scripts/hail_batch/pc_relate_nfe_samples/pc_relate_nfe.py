"""
Perform pc_relate on nfe samples from the HGDP/1KG dataset.
"""

import hail as hl
from analysis_runner import output_path
from hail.experimental import lgt_to_gt

TOB_WGS = 'gs://cpg-tob-wgs-test/densify/v1/tob_wgs_filtered.mt/'


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(TOB_WGS)
    mt = mt.annotate_entries(GT=lgt_to_gt(mt.LGT, mt.LA))
    nrows_mt = mt.count_rows()
    mt = mt.sample_rows(100 / nrows_mt, seed=12345)

    # Remove related samples (at the 2nd degree or closer)
    king = hl.king(mt.GT)
    known_trios = king.filter_cols((king.s == 'HG01696') | (king.s == 'HG01629'))
    # Select parents from the child_pca_outlier matrix
    known_trios = known_trios.filter_rows(
        (known_trios.s_1 == 'HG01628')
        | (known_trios.s_1 == 'HG01694')
        | (known_trios.s_1 == 'HG01695')
        | (known_trios.s_1 == 'HG01630')
    )
    print(known_trios.entries().to_pandas())
    king_path = output_path('king_kinship_estimate.ht')
    known_trios.entries().write(king_path)


if __name__ == '__main__':
    query()
