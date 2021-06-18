"""
Perform pca on samples specific to a population
from the HGDP,1KG, and tob-wgs dataset after densifying.

Depends on hgdp1kg_tobwgs_densified_pca/hgdp_1kg_tob_wgs_densified_pca.py
"""

import click
import pandas as pd
import hail as hl

HGDP1KG_TOBWGS = (
    'gs://cpg-tob-wgs-main/1kg_hgdp_densified_pca/v2/'
    'hgdp1kg_tobwgs_joined_all_samples.mt'
)


@click.command()
@click.option('--output', help='GCS output path', required=True)
@click.option('--pop', help='Population to subset from the 1KG (e.g. afr, nfe)')
def query(output, pop):  # pylint: disable=too-many-locals
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(HGDP1KG_TOBWGS)
    if pop:
        # Get samples from the specified population only
        mt = mt.filter_cols(
            (mt.hgdp_1kg_metadata.population_inference.pop == pop.lower())
            | (mt.s.contains('TOB'))
        )
    else:
        mt = mt.filter_cols(mt.s.contains('TOB'))

    # Perform PCA
    eigenvalues_path = f'{output}/eigenvalues.ht'
    scores_path = f'{output}/scores.ht'
    loadings_path = f'{output}/loadings.ht'
    eigenvalues, scores, loadings = hl.hwe_normalized_pca(
        mt.GT, compute_loadings=True, k=20
    )
    hl.Table.from_pandas(pd.DataFrame(eigenvalues)).export(eigenvalues_path)
    scores.write(scores_path, overwrite=True)
    loadings.write(loadings_path, overwrite=True)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
