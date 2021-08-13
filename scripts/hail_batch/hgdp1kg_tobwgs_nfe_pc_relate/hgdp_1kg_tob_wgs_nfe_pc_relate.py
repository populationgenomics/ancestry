"""
Estimate kinship for nfe samples from the HGDP/1KG +
tob-wgs datasets.

"""

import hail as hl
import pandas as pd
from analysis_runner import bucket_path, output_path


HGDP1KG_TOBWGS = bucket_path(
    '1kg_hgdp_densified_pca_new_variants/v0/hgdp1kg_tobwgs_joined_all_samples.mt'
)


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(HGDP1KG_TOBWGS)
    # Get samples from the specified population only
    mt = mt.filter_cols(
        (mt.hgdp_1kg_metadata.population_inference.pop == 'nfe')
        | (mt.s.contains('TOB'))
    )

    # Perform kinship test with pc_relate
    pc_rel = hl.pc_relate(mt.GT, 0.01, k=20, statistics='kin')
    child_pca_outliers = pc_rel.filter(
        (pc_rel.i.s == 'HG01696') | (pc_rel.i.s == 'HG01629')
    )
    # Select parents from the child_pca_outlier matrix
    child_pca_outliers = child_pca_outliers.filter(
        (child_pca_outliers.j.s == 'HG01628')
        | (child_pca_outliers.j.s == 'HG01694')
        | (child_pca_outliers.j.s == 'HG01695')
        | (child_pca_outliers.j.s == 'HG01630')
    )

    # save as html
    html = pd.DataFrame(
        {
            'individual_i': child_pca_outliers.i.s.collect(),
            'individual_j': child_pca_outliers.j.s.collect(),
            # this number ranges from 0 to 0.5, with 0.5 being identical
            'kinship_estimate': child_pca_outliers.kin.collect(),
        }
    ).to_html()
    plot_filename_html = output_path(f'sample_outliers_pc_relate_estimates.html', 'web')
    with hl.hadoop_open(plot_filename_html, 'w') as f:
        f.write(html)


if __name__ == '__main__':
    query()
