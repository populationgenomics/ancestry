"""Test cost and speed of densifying on 20 samples"""

import hail as hl
from analysis_runner import bucket_path, output_path

TOB_WGS = bucket_path('mt/v5.1.mt')


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(TOB_WGS)
    # save output
    mt_path = output_path(f'sparse_tob_wgs.mt', 'tmp')
    mt.write(mt_path)


if __name__ == '__main__':
    query()
