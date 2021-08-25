"""
Estimate kinship for samples from the HGDP/1KG +
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

    # Perform kinship test with pc_relate
    pc_rel_path = output_path('pc_relate_kinship_estimate.ht')
    pc_rel = hl.pc_relate(mt.GT, 0.01, k=10, statistics='kin')
    pc_rel.write(pc_rel_path, overwrite=True)
    pairs = pc_rel.filter(pc_rel['kin'] >= 0.125)
    related_samples_to_remove = hl.maximal_independent_set(pairs.i, pairs.j, False)
    n_related_samples = related_samples_to_remove.count()
    print(f'related_samples_to_remove.count() = {n_related_samples}')

    # save as html
    html = pd.DataFrame(
        {'removed_individual': related_samples_to_remove.node.s.collect()}
    ).to_html()
    plot_filename_html = output_path(f'removed_samples.html', 'web')
    with hl.hadoop_open(plot_filename_html, 'w') as f:
        f.write(html)


if __name__ == '__main__':
    query()
