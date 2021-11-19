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
    'gs://cpg-tob-wgs-main-analysis/joint-calling/v7/ancestry/mt_union_hgdp.mt'
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
        'HG01628',
        'NA12892',
        'HG01694',
        'HG01695',
        'HG01669',
        'HG01527',
        'HG01630',
        'CPG10066',
        'CPG1412',
        'CPG1677',
        'CPG1784',
        'CPG1842',
        'CPG2329',
        'CPG2691',
        'CPG3129',
        'CPG3350',
        'CPG3467',
        'CPG4374',
        'CPG6122',
        'CPG8250',
        'CPG9191',
        'CPG9704',
        'HG01630',
        'CPG8532',
        'CPG1099',
        'CPG3822',
        'CPG2519',
        'HG01697',
        'HG01757',
        'HG01704',
        'v3.1::HG02232',
        'HG01668',
        'HG01612',
        'HG01766',
        'HG01506',
        'CPG7633',
        'CPG1248',
        'CPG1347',
        'CPG3731',
        'HG01626',
        'HG01500',
        'HG01619',
        'HGDP01384',
        'HG01767',
        'HG01525',
        'CPG919',
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
