"""
Generate PCA on combined TOB-WGS/SNP-chip data.
"""

import click
import hail as hl
import pandas as pd

TOB_WGS = 'gs://cpg-tob-wgs-main/mt/v2-raw.mt/'
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
    tob_wgs = tob_wgs.select_cols().select_entries(tob_wgs.GT)
    snp_chip = snp_chip.select_cols().select_entries(snp_chip.GT)
    tob_combined = tob_wgs.union_cols(snp_chip)
    tob_combined = tob_combined.cache()
    print(tob_combined.count_rows())
    tob_combined_path = f'{output}/tob_wgs_snp_chip_combined.mt'
    tob_combined = tob_wgs.repartition(1000, shuffle=False)
    tob_combined.write(tob_combined_path)

    # Perform PCA
    eigenvalues_path = f'{output}/eigenvalues.ht'
    scores_path = f'{output}/scores.ht'
    loadings_path = f'{output}/loadings.ht'
    eigenvalues, scores, loadings = hl.hwe_normalized_pca(
        tob_combined.GT, compute_loadings=True, k=20
    )
    hl.Table.from_pandas(pd.DataFrame(eigenvalues)).export(eigenvalues_path)
    scores.write(scores_path, overwrite=True)
    loadings.write(loadings_path, overwrite=True)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
