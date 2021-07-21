"""
Increase the number of partitions
in the SNP-Chip dataset.
"""

import hail as hl
from analysis_runner import bucket_path, output_path

SNP_CHIP = bucket_path('snpchip/v1/snpchip_grch38.mt')
# TOB_WGS = bucket_path('mt/v3-raw.mt')
TOB_WGS = 'gs://cpg-tob-wgs-test/mt/v4.mt/'


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    snp_chip = hl.read_matrix_table(SNP_CHIP).key_rows_by('locus', 'alleles')
    tob_wgs = hl.read_matrix_table(TOB_WGS)
    snp_chip = snp_chip.semi_join_rows(tob_wgs.rows())
    snp_chip = snp_chip.repartition(100)
    snp_chip_path = output_path('snp_chip_100_partitions.mt')
    snp_chip.write(snp_chip_path, overwrite=True)


if __name__ == '__main__':
    query()
