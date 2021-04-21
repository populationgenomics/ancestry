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

    mt = hl.read_matrix_table(GNOMAD_HGDP_1KG_MT)
    eigenvalues_path = f'{output}/eigenvalues_10k.csv'
    scores_path = f'{output}/scores_10k.ht'
    loadings_path = f'{output}/loadings_10k.ht'

    # filter out variants with a call rate <0.99
    mt_qc = hl.variant_qc(mt)
    filt_mt = mt_qc.filter_rows(mt_qc.variant_qc.call_rate >= 0.99)
    nrows = filt_mt.count_rows()
    # Downsample the dataset to approximately 10k randomly-selected rows
    # (the input must be a proportion)
    downsampled_mt = filt_mt.sample_rows(10000 / nrows)

    # test on 100 samples
    mt_head = downsampled_mt.head(None, n_cols=100)
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
