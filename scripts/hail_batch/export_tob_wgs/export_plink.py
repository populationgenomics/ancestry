"""Export TOB-WGS joint callset as PLINK format"""

import hail as hl
from analysis_runner import bucket_path, output_path

TOB_WGS = bucket_path('mt/v5.1.mt/')


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    tob_wgs = hl.read_matrix_table(TOB_WGS)
    tob_wgs = hl.experimental.densify(tob_wgs)
    tob_wgs = hl.split_multi_hts(tob_wgs)
    tob_wgs_path = output_path('tob_wgs_plink')
    hl.export_plink(tob_wgs, tob_wgs_path, ind_id=tob_wgs.s)


if __name__ == '__main__':
    query()
