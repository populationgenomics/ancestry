"""
Perform pc_relate on nfe samples from the HGDP/1KG dataset.
"""

import hail as hl
import pandas as pd
from analysis_runner import bucket_path, output_path

KING_GLOBAL = bucket_path('pc_relate/v3/king_kinship_estimate_global.ht')
KING_NFE = bucket_path('pc_relate/v2/king_kinship_estimate.ht')


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    # Remove related samples (at the 2nd degree or closer)
    king = hl.read_matrix_table(KING_GLOBAL)
    ht = king.entries()
    related_samples = ht.filter((ht.s_1 != ht.s) & (ht.phi > 0.1), keep=True)
    king_global = pd.DataFrame(
        {
            'i_s': related_samples.s_1.collect(),
            'j_s': related_samples.s.collect(),
            'kin': related_samples.phi.collect(),
        }
    )
    filename = output_path('king_global_matrix_10k.csv')
    king_global.to_csv(filename, index=False)

    # get maximal independent set
    second_degree_related_samples = ht.filter(
        (ht.s_1 != ht.s) & (ht.phi > 0.125), keep=True
    )
    struct = hl.struct(
        i=second_degree_related_samples.s_1, j=second_degree_related_samples.s
    )
    struct = struct.annotate(phi=second_degree_related_samples.phi)
    related_samples_to_remove = hl.maximal_independent_set(
        struct.i, struct.j, False  # pylint: disable=E1101
    )
    related_samples = pd.DataFrame(
        {'related_individual': related_samples_to_remove.node.collect()}
    )
    filename = output_path(
        'king_10k_global_related_samples_maximal_independent_set.csv'
    )
    related_samples.to_csv(filename, index=False)

    # Remove related samples (at the 2nd degree or closer) for NFE samples
    king_nfe = hl.read_matrix_table(KING_NFE)
    ht = king_nfe.entries()
    related_samples = ht.filter((ht.s_1 != ht.s) & (ht.phi > 0.1), keep=True)
    king_nfe = pd.DataFrame(
        {
            'i_s': related_samples.s_1.collect(),
            'j_s': related_samples.s.collect(),
            'kin': related_samples.phi.collect(),
        }
    )
    filename = output_path('king_nfe_matrix_10k.csv')
    king_nfe.to_csv(filename, index=False)

    # get maximal independent set
    second_degree_related_samples = ht.filter(
        (ht.s_1 != ht.s) & (ht.phi > 0.125), keep=True
    )
    struct = hl.struct(
        i=second_degree_related_samples.s_1, j=second_degree_related_samples.s
    )
    struct = struct.annotate(phi=second_degree_related_samples.phi)
    related_samples_to_remove = hl.maximal_independent_set(
        struct.i, struct.j, False  # pylint: disable=E1101
    )
    related_samples = pd.DataFrame(
        {'related_individual': related_samples_to_remove.node.collect()}
    )
    filename = output_path('king_10k_nfe_related_samples_maximal_independent_set.csv')
    related_samples.to_csv(filename, index=False)


if __name__ == '__main__':
    query()
