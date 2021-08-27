"""
Estimate kinship coefficient using KING on TOB samples from the HGDP/1KG dataset.
"""

import hail as hl
from analysis_runner import bucket_path, output_path

HGDP1KG_TOBWGS = bucket_path(
    '1kg_hgdp_densified_pca_new_variants/v0/hgdp1kg_tobwgs_joined_all_samples.mt'
)


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(HGDP1KG_TOBWGS)
    mt = mt.filter_cols(mt.s.contains('TOB'))
    # Remove related samples (at the 2nd degree or closer)
    king = hl.king(mt.GT)
    king_path = output_path('king_kinship_estimate_TOB.ht')
    king.write(king_path)


if __name__ == '__main__':
    query()
