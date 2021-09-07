"""
Save scores of related individuals after running pc_relate.

"""

import hail as hl
import pandas as pd
from analysis_runner import bucket_path, output_path


PC_RELATE_ESTIMATE_NFE = bucket_path(
    'tob_wgs_hgdp_1kg_nfe_pc_relate/v0/pc_relate_kinship_estimate.ht'
)
PC_RELATE_ESTIMATE_GLOBAL = bucket_path(
    'tob_wgs_hgdp_1kg_pc_relate/v0/pc_relate_kinship_estimate.ht'
)
KING_ESTIMATE_NFE = bucket_path('king/v0/king_kinship_estimate_NFE.ht')


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    # save relatedness estimates for pc_relate global populations
    ht = hl.read_table(PC_RELATE_ESTIMATE_GLOBAL)
    related_samples = ht.filter(ht.kin > 0.1)
    pc_relate_global = pd.DataFrame(
        {
            'i_s': related_samples.i.s.collect(),
            'j_s': related_samples.j.s.collect(),
            'kin': related_samples.kin.collect(),
        }
    )
    filename = output_path(f'pc_relate_global_matrix.csv', 'analysis')
    pc_relate_global.to_csv(filename, index=False)

    # get maximal independent set
    pairs = ht.filter(ht['kin'] >= 0.125)
    related_samples_to_remove = hl.maximal_independent_set(pairs.i, pairs.j, False)

    related_samples = pd.DataFrame(
        {'removed_individual': related_samples_to_remove.node.s.collect()}
    )
    filename = output_path(f'pc_relate_global_maximal_independent_set.csv', 'analysis')
    related_samples.to_csv(filename, index=False)

    # save relatedness estimates for pc_relate NFE samples
    ht = hl.read_table(PC_RELATE_ESTIMATE_NFE)
    related_samples = ht.filter(ht.kin > 0.1)
    pc_relate_nfe = pd.DataFrame(
        {
            'i_s': related_samples.i.s.collect(),
            'j_s': related_samples.j.s.collect(),
            'kin': related_samples.kin.collect(),
        }
    )
    filename = output_path(f'pc_relate_nfe_matrix.csv', 'analysis')
    pc_relate_nfe.to_csv(filename, index=False)
    # get maximal independent set
    pairs = ht.filter(ht['kin'] >= 0.125)
    related_samples_to_remove = hl.maximal_independent_set(pairs.i, pairs.j, False)
    related_samples = pd.DataFrame(
        {'removed_individual': related_samples_to_remove.node.s.collect()}
    )
    filename = output_path(f'pc_relate_nfe_maximal_independent_set.csv', 'analysis')
    related_samples.to_csv(filename, index=False)

    # save relatedness estimates for KING NFE samples
    mt = hl.read_matrix_table(KING_ESTIMATE_NFE)
    ht = mt.entries()
    # remove entries where samples are identical
    related_samples = ht.filter(ht.s_1 != ht.s)
    related_samples = ht.filter(ht.phi > 0.1)
    king_nfe = pd.DataFrame(
        {
            'i_s': related_samples.s_1.collect(),
            'j_s': related_samples.s.collect(),
            'kin': related_samples.phi.collect(),
        }
    )
    filename = output_path(f'king_nfe_matrix_90k.csv', 'analysis')
    king_nfe.to_csv(filename, index=False)
    # save KING NFE maximal independent set
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
        f'king_90k_related_samples_maximal_independent_set.csv', 'analysis'
    )
    related_samples.to_csv(filename, index=False)


if __name__ == '__main__':
    query()
