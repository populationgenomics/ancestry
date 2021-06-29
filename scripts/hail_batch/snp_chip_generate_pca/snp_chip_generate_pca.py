"""
Test densify function on TOB-WGS data.
"""

import click
import hail as hl
import pandas as pd

# TOB_WGS = 'gs://cpg-tob-wgs-main/mt/v2-raw.mt/'
TOB_WGS = 'gs://cpg-tob-wgs-test/mt/v3.mt/'
SNP_CHIP = 'gs://cpg-tob-wgs-test/snpchip/v1/snpchip_grch38.mt/'


@click.command()
@click.option('--output', help='GCS output path', required=True)
def query(output):  # pylint: disable=too-many-locals
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    tob_wgs = hl.read_matrix_table(TOB_WGS).key_rows_by('locus', 'alleles')
    snp_chip = hl.read_matrix_table(SNP_CHIP).key_rows_by('locus', 'alleles')

    # filter to loci that are contained in snp-chip data after densifying
    tob_wgs = hl.experimental.densify(tob_wgs)
    tob_wgs = tob_wgs.head(100000)
    tob_wgs = tob_wgs.semi_join_rows(snp_chip.rows())
    tob_wgs = tob_wgs.cache()
    print(tob_wgs.count_rows())
    tob_wgs_path = f'{output}/tob_wgs_filtered_snp_chip.mt'
    tob_wgs = tob_wgs.repartition(100, shuffle=False)
    tob_wgs.write(tob_wgs_path)

    # Perform PCA
    eigenvalues_path = f'{output}/eigenvalues.ht'
    scores_path = f'{output}/scores.ht'
    loadings_path = f'{output}/loadings.ht'
    eigenvalues, scores, loadings = hl.hwe_normalized_pca(
        tob_wgs.GT, compute_loadings=True, k=20
    )
    hl.Table.from_pandas(pd.DataFrame(eigenvalues)).export(eigenvalues_path)
    scores.write(scores_path, overwrite=True)
    loadings.write(loadings_path, overwrite=True)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
