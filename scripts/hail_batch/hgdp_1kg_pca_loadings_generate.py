"""Generates loadings, scores, and eigenvalues for the HGDP + 1KG dataset"""

import click
import pandas as pd
import hail as hl

GNOMAD_HGDP_1KG_MT = (
    'gs://gcp-public-data--gnomad/release/3.1/mt/genomes/'
    'gnomad.genomes.v3.1.hgdp_1kg_subset_dense.mt'
)


@click.command()
@click.option('--output', help='GCS output path', required=True)
def query(output):
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    eigenvalues_path = f'{output}/eigenvalues.csv'
    scores_path = f'{output}/scores.ht'
    loadings_path = f'{output}/loadings.ht'
    mt = hl.read_matrix_table(GNOMAD_HGDP_1KG_MT)
    # test on 100 samples
    mt_head = mt.head(n=mt.count_rows(), n_cols=20)
    eigenvalues, scores, loadings = hl.hwe_normalized_pca(
        mt_head.GT, compute_loadings=True, k=20
    )
    # save the list of eigenvalues
    eigenvalues_df = pd.DataFrame(eigenvalues)
    eigenvalues_df.to_csv(eigenvalues_path, index=False)
    # save the scores and loadings as a hail table
    scores.write(scores_path, overwrite=True)
    loadings.write(loadings_path, overwrite=True)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
