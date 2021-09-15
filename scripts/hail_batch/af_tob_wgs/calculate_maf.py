"""Get number of variants with a MAF < 0.05"""

import hail as hl
from analysis_runner import bucket_path, output_path

TOB_WGS = bucket_path('mt/v5.1.mt/')


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    tob_wgs = hl.read_matrix_table(TOB_WGS)
    tob_wgs = hl.experimental.densify(tob_wgs)
    tob_wgs = hl.variant_qc(tob_wgs)
    # get MAF < 0.05
    snp_maf_05 = tob_wgs.aggregate_rows(
        hl.agg.count_where(tob_wgs.variant_qc.AF[1] < 0.05)
    )
    print(f'Variant MAF < 0.05 = {snp_maf_05}')


if __name__ == '__main__':
    query()
