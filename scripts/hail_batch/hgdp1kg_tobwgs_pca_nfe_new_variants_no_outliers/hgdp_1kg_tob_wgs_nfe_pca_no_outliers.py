"""
Perform pca on nfe samples from the HGDP/1KG +
tob-wgs dataset and remove outlier samples.
Reliant on output from

```hgdp1kg_tobwgs_densified_pca_new_variants/
hgdp_1kg_tob_wgs_densified_pca_new_variants.py
```
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
    # remove outlier samples, as identified by PCA
    outliers = [
        'TOB1734',
        'TOB1714',
        'TOB1126',
        'TOB1653',
        'TOB1668',
        'TOB1681',
        'TOB1116',
        'TOB1107',
        'TOB1635',
        'HG01628',
        'TOB1675',
        'TOB1125',
        'TOB1762',
        'TOB1263',
        'TOB1640',
        'HG01669',
        'TOB1795',
        'TOB1707',
        'HG01695',
        'HG01694',
        'TOB1673',
        'HG01630',
    ]

    mt = mt.filter_cols(hl.literal(outliers).contains(mt.s), keep=False)

    # Remove related samples at the 2nd degree or closer, as indicated by gnomAD
    mt = mt.filter_cols(mt.hgdp_1kg_metadata.gnomad_release | mt.s.startswith('TOB'))

    # Perform PCA
    eigenvalues_path = output_path('eigenvalues.ht')
    scores_path = output_path('scores.ht')
    loadings_path = output_path('loadings.ht')
    eigenvalues, scores, loadings = hl.hwe_normalized_pca(
        mt.GT, compute_loadings=True, k=20
    )
    hl.Table.from_pandas(pd.DataFrame(eigenvalues)).export(eigenvalues_path)
    scores.write(scores_path, overwrite=True)
    loadings.write(loadings_path, overwrite=True)


if __name__ == '__main__':
    query()
