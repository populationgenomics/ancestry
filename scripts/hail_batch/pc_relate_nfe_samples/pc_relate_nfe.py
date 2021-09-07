"""
Perform pc_relate on nfe samples from the HGDP/1KG dataset.
"""

import hail as hl
import pandas as pd
from analysis_runner import bucket_path, output_path

KING = bucket_path('pc_relate/v3/king_kinship_estimate_global.ht')


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    # Remove related samples (at the 2nd degree or closer)
    king = hl.read_matrix_table(KING)
    ht = king.entries()
    related_samples = ht.filter((ht.s_1 != ht.s) & (ht.phi > 0.1), keep=True)
    king_global = pd.DataFrame(
        {
            'i_s': related_samples.s_1.collect(),
            'j_s': related_samples.s.collect(),
            'kin': related_samples.phi.collect(),
        }
    )
    filename = output_path(f'king_global_matrix_90k.csv', 'metadata')
    king_global.to_csv(filename, index=False)


if __name__ == '__main__':
    query()
