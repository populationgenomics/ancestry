"""
Perform pca on nfe samples from the HGDP/1KG +
tob-wgs dataset and remove outlier samples.
Reliant on output from
```joint-calling/scripts/ancestry_pca.py```
"""

import hail as hl
import pandas as pd
from analysis_runner import output_path


HGDP1KG_TOBWGS = (
    'gs://cpg-tob-wgs-main-analysis/joint-calling/v6/ancestry/mt_union_hgdp.mt'
)


def query():
    """Query script entry point."""
    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(HGDP1KG_TOBWGS)
    # Get samples from the specified population only
    mt = mt.filter_cols(
        (mt.hgdp_1kg_metadata.population_inference.pop == 'nfe')
        | (mt.s.contains('CPG'))
    )
    # remove outlier samples, as identified by PCA
    outliers = [
        'CPG2709',
        'CPG2881',
        'CPG3707',
        'CPG1040',
        'CPG2261',
        'CPG9175',
        'CPG9142',
        'CPG4804',
        'CPG2147',
        'CPG3517',
        'CPG3608',
        'CPG2386',
        'CPG2337',
        'CPG3699',
        'CPG786',
        'CPG927',
        'CPG3244',
        'CPG2634',
        'CPG2311',
        'CPG6122',
        'CPG1677',
        'CPG9191',
        'CPG9704',
        'CPG2329',
        'CPG2691',
        'CPG8250',
        'CPG3129',
        'CPG3467',
        'CPG3350',
        'CPG10066',
        'CPG1784',
        'CPG4374',
        'CPG1412',
        'CPG1099',
        'CPG1842',
    ]

    mt = mt.filter_cols(hl.literal(outliers).contains(mt.s), keep=False)

    # Remove related samples at the 2nd degree or closer, as indicated by gnomAD
    mt = mt.filter_cols(mt.hgdp_1kg_metadata.gnomad_release | mt.s.startswith('CPG'))

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
