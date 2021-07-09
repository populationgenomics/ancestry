"""
Generate PCA on SNP-chip data only.
"""

import click
import hail as hl
import pandas as pd
from analysis_runner import bucket_path, output_path

SNP_CHIP = bucket_path('snpchip/v1/snpchip_grch38.mt')


@click.command()
def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    snp_chip = hl.read_matrix_table(SNP_CHIP)
    eigenvalues_path = output_path('eigenvalues.ht')
    scores_path = output_path('scores.ht')
    loadings_path = output_path('loadings.ht')
    # Perform PCA
    eigenvalues, scores, loadings = hl.hwe_normalized_pca(
        snp_chip.GT, compute_loadings=True, k=5
    )
    hl.Table.from_pandas(pd.DataFrame(eigenvalues)).export(eigenvalues_path)
    scores.write(scores_path, overwrite=True)
    loadings.write(loadings_path, overwrite=True)


if __name__ == '__main__':
    query()
