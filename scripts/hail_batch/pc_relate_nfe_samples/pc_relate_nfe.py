"""
Perform pc_relate on nfe samples from the HGDP/1KG dataset.
"""

import hail as hl
from analysis_runner import output_path

GNOMAD_HGDP_1KG_MT = (
    'gs://gcp-public-data--gnomad/release/3.1/mt/genomes/'
    'gnomad.genomes.v3.1.hgdp_1kg_subset_dense.mt'
)


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(GNOMAD_HGDP_1KG_MT)
    nrows_mt = mt.count_rows()
    mt = mt.sample_rows(10000 / nrows_mt, seed=12345)
    mt_path = output_path('hgdp_1kg_subsampled_10k.mt', 'tmp')
    mt = mt.checkpoint(mt_path)
    mt = mt.repartition(1000, shuffle=False)
    mt = mt.filter_cols(mt.population_inference.pop == 'nfe')
    mt_repartitioned_path = output_path('hgdp_1kg_repartitioned_10k.mt', 'tmp')
    mt = mt.checkpoint(mt_repartitioned_path)

    # Remove related samples (at the 2nd degree or closer)
    king = hl.king(mt.GT)
    king_path = output_path('king_kinship_estimate.ht')
    king.write(king_path)


if __name__ == '__main__':
    query()
