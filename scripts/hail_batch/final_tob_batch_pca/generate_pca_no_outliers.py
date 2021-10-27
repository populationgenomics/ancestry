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
    'gs://cpg-tob-wgs-main-analysis/joint-calling/v6-1/ancestry/mt_union_hgdp.mt'
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
        'CPG9142',
        'CPG9191',
        'CPG9704',
        'HG01694',
        'HG01695',
        'HG01669',
        'HG01527',
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
        #     'CPG1206',
        #     'CPG4168',
        #     'CPG5611',
        #     'CPG1925',
        #     'CPG7328',
        #     'CPG2550',
        #     'CPG869',
        #     'CPG10355',
        #     'CPG10603',
        #     'CPG1149',
        #     'CPG1644',
        #     'CPG1701',
        #     'CPG1727',
        #     'CPG1875',
        #     'CPG2014',
        #     'CPG2071',
        #     'CPG3004',
        #     'CPG3665',
        #     'CPG737',
        #     'CPG1420',
        #     'CPG1537',
        #     'CPG1560',
        #     'CPG1891',
        #     'CPG1974',
        #     'CPG2022',
        #     'CPG2055',
        #     'CPG2402',
        #     'CPG3186',
        #     'CPG3236',
        #     'CPG3574',
        #     'CPG3590',
        #     'CPG7120',
        #     'CPG1032',
        #     'CPG1198',
        #     'CPG1271',
        #     'CPG1362',
        #     'CPG1776',
        #     'CPG1982',
        #     'CPG2378',
        #     'CPG2444',
        #     'CPG2840',
        #     'CPG2907',
        #     'CPG3418',
        #     'CPG3756',
        #     'CPG3772',
        #     'CPG4176',
        #     'CPG794',
        #     'CPG810',
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
