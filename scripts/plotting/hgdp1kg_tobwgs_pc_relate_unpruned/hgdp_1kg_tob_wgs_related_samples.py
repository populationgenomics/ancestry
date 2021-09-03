"""
Save scores of related individuals after running pc_relate.

"""

import hail as hl
import pandas as pd
from analysis_runner import bucket_path, output_path


# KINSHIP_ESTIMATE_NFE = bucket_path(
#     'tob_wgs_hgdp_1kg_nfe_pc_relate/v0/pc_relate_kinship_estimate.ht'
# )

KINSHIP_ESTIMATE_NFE = bucket_path(
    'tob_wgs_hgdp_1kg_nfe_pc_relate/v1/pc_relate_kinship_estimate.ht'
)

# KINSHIP_ESTIMATE_GLOBAL = bucket_path(
#     'tob_wgs_hgdp_1kg_pc_relate/v0/pc_relate_kinship_estimate.ht'
# )

KINSHIP_ESTIMATE_GLOBAL = bucket_path(
    'tob_wgs_hgdp_1kg_nfe_pc_relate/v1/pc_relate_kinship_estimate.ht'
)

# KING_ESTIMATE_NFE = bucket_path('king/v0/king_kinship_estimate_NFE.ht')

KING_ESTIMATE_NFE = bucket_path('pc_relate/v3/king_kinship_estimate_global.ht')


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    # save relatedness estimates for pc_relate global populations
    ht = hl.read_table(KINSHIP_ESTIMATE_GLOBAL)
    pc_relate_global = pd.DataFrame(
        {
            'i_s': ht.i.s.collect(),
            'j_s': ht.j.s.collect(),
            'kin': ht.kin.collect(),
        }
    )
    filename = output_path(f'pc_relate_global_matrix.csv', 'metadata')
    pc_relate_global.to_csv(filename, index=False)

    # save relatedness estimates for pc_relate NFE samples
    ht = hl.read_table(KINSHIP_ESTIMATE_NFE)
    pc_relate_nfe = pd.DataFrame(
        {
            'i_s': ht.i.s.collect(),
            'j_s': ht.j.s.collect(),
            'kin': ht.kin.collect(),
        }
    )
    filename = output_path(f'pc_relate_nfe_matrix.csv', 'metadata')
    pc_relate_nfe.to_csv(filename, index=False)

    # save relatedness estimates for KING NFE samples
    mt = hl.read_matrix_table(KING_ESTIMATE_NFE)
    ht = mt.entries()
    # remove entries where samples are identical
    related_samples = ht.filter(ht.s_1 != ht.s)
    king_nfe = pd.DataFrame(
        {
            'i_s': related_samples.s_1.collect(),
            'j_s': related_samples.s.collect(),
            'kin': related_samples.phi.collect(),
        }
    )
    filename = output_path(f'king_nfe_matrix_90k.csv', 'metadata')
    king_nfe.to_csv(filename, index=False)


if __name__ == '__main__':
    query()
