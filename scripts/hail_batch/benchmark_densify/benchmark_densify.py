"""Test cost and speed of densifying on 20 samples"""

import hail as hl
from analysis_runner import bucket_path, output_path

TOB_WGS = bucket_path('mt/v5.1.mt')


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(TOB_WGS)
    # remove unecessary information from the mt
    densified_mt = hl.experimental.densify(mt)
    mt = mt.select_entries(mt.GT).select_globals().select_rows().select_cols()
    densified_mt = (
        densified_mt.select_entries(densified_mt.GT)
        .select_globals()
        .select_rows()
        .select_cols()
    )
    # save output
    mt_path = output_path(f'sparse_tob_wgs.mt', 'tmp')
    densified_path = output_path(f'densified_tob_wgs.mt', 'tmp')
    mt.write(mt_path)
    densified_mt.write(densified_path)


if __name__ == '__main__':
    query()
