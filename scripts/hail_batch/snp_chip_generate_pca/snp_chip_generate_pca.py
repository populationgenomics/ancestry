"""
Generate PCA on combined TOB-WGS/SNP-chip data.
"""

import click
import hail as hl
import pandas as pd
from hail.experimental import lgt_to_gt
from analysis_runner import output_path

TOB_WGS = 'gs://cpg-tob-wgs-main/mt/v3-raw.mt/'
SNP_CHIP = 'gs://cpg-tob-wgs-test/snpchip/v1/snpchip_grch38.mt/'


@click.command()
def query():  # pylint: disable=too-many-locals
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    tob_wgs = hl.read_matrix_table(TOB_WGS).key_rows_by('locus', 'alleles')
    snp_chip = hl.read_matrix_table(SNP_CHIP).key_rows_by('locus', 'alleles')

    # filter to loci that are contained in snp-chip data after densifying
    tob_wgs = hl.experimental.densify(tob_wgs)
    tob_wgs = tob_wgs.select_entries(
        GT=lgt_to_gt(tob_wgs.LGT, tob_wgs.LA)
    ).select_cols()
    snp_chip = snp_chip.select_entries(snp_chip.GT).select_cols()
    snp_chip = snp_chip.key_cols_by(s=snp_chip.s + '_snp_chip')
    tob_combined = tob_wgs.union_cols(snp_chip)
    tob_combined = tob_combined.cache()
    print(tob_combined.count_rows())

    # Perform PCA
    eigenvalues_path = output_path('eigenvalues.ht')
    scores_path = output_path('scores.ht')
    loadings_path = output_path('loadings.ht')
    eigenvalues, scores, loadings = hl.hwe_normalized_pca(
        tob_combined.GT, compute_loadings=True, k=20
    )
    hl.Table.from_pandas(pd.DataFrame(eigenvalues)).export(eigenvalues_path)
    scores.write(scores_path, overwrite=True)
    loadings.write(loadings_path, overwrite=True)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
