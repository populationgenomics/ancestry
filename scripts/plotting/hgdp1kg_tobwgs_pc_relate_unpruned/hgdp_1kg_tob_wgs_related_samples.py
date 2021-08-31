"""
Plot scores of related individuals after running pc_relate.

"""

import hail as hl
import pandas as pd
from analysis_runner import bucket_path, output_path


KINSHIP_ESTIMATE_NFE = bucket_path(
    'tob_wgs_hgdp_1kg_nfe_pc_relate/v0/pc_relate_kinship_estimate.ht'
)

KINSHIP_ESTIMATE_GLOBAL = bucket_path(
    'tob_wgs_hgdp_1kg_pc_relate/v0/pc_relate_kinship_estimate.ht'
)


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    # plot relatedness estimates for global populations
    ht = hl.read_table(KINSHIP_ESTIMATE_GLOBAL)

    related_samples = ht.filter(ht.kin > 0.1)

    # save as html
    html = pd.DataFrame(
        {
            'i_s': related_samples.i.s.collect(),
            'j_s': related_samples.j.s.collect(),
            'kin': related_samples.kin.collect(),
        }
    ).to_html()
    plot_filename_html = output_path(f'pc_relate_global_matrix.html', 'web')
    with hl.hadoop_open(plot_filename_html, 'w') as f:
        f.write(html)

    # plot relatedness estimates for NFE samples
    ht = hl.read_table(KINSHIP_ESTIMATE_NFE)

    related_samples = ht.filter(ht.kin > 0.1)

    # save as html
    html = pd.DataFrame(
        {
            'i_s': related_samples.i.s.collect(),
            'j_s': related_samples.j.s.collect(),
            'kin': related_samples.kin.collect(),
        }
    ).to_html()
    plot_filename_html = output_path(f'pc_relate_nfe_matrix.html', 'web')
    with hl.hadoop_open(plot_filename_html, 'w') as f:
        f.write(html)


if __name__ == '__main__':
    query()
