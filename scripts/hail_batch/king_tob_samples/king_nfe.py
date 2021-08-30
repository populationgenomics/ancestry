"""
Estimate kinship coefficient using KING on NFE samples from the HGDP/1KG dataset.
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
    mt = mt.filter_cols(
        (mt.hgdp_1kg_metadata.population_inference.pop == 'nfe')
        | (mt.s.contains('TOB'))
    )
    # Remove related samples (at the 2nd degree or closer)
    king = hl.king(mt.GT)
    king_path = output_path('king_kinship_estimate_NFE.ht')
    king.write(king_path)
    ht = king.entries()
    related_samples = ht.filter((ht.s_1 != ht.s) & (ht.phi > 0.125), keep=True)
    struct = hl.struct(i=related_samples.s_1, j=related_samples.s)
    struct = struct.annotate(phi=related_samples.phi)
    related_samples_to_remove = hl.maximal_independent_set(
        struct.i, struct.j, False  # pylint: disable=E1101
    )
    n_related_samples = related_samples_to_remove.count()
    print(f'related_samples_to_remove.count() = {n_related_samples}')
    # save as html
    html = pd.DataFrame(
        {'related_individual': related_samples_to_remove.node.collect()}
    ).to_html()
    plot_filename_html = output_path(f'related_samples.html', 'web')
    with hl.hadoop_open(plot_filename_html, 'w') as f:
        f.write(html)


if __name__ == '__main__':
    query()
