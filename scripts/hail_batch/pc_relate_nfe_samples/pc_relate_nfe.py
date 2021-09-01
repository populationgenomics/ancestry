"""
Perform pc_relate on nfe samples from the HGDP/1KG dataset.
"""

import hail as hl
from analysis_runner import output_path

HGDP_1KG_MT_SUBSAMPLED = (
    'gs://cpg-tob-wgs-test-tmp/pc_relate/v2/hgdp_1kg_subsampled_10k.mt/'
)


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(HGDP_1KG_MT_SUBSAMPLED)
    mt = mt.repartition(1000, shuffle=False)
    mt_repartitioned_path = output_path('hgdp_1kg_repartitioned_10k_global.mt', 'tmp')
    mt = mt.checkpoint(mt_repartitioned_path)

    # Remove related samples (at the 2nd degree or closer)
    king = hl.king(mt.GT)
    king_path = output_path('king_kinship_estimate_global.ht')
    king.write(king_path)


if __name__ == '__main__':
    query()
