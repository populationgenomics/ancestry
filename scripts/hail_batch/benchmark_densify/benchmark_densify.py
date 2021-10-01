"""Test cost and speed of densifying on 20 samples"""

import hail as hl
from analysis_runner import output_path

TOB_WGS = 'gs://cpg-tob-wgs-test-tmp/joint-calling/v6-25/combiner/v6-25-raw.mt'


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(TOB_WGS)
    mt = hl.experimental.densify(mt)
    # save output
    mt_path = output_path(f'densified_tob_wgs.mt', 'tmp')
    mt.write(mt_path)


if __name__ == '__main__':
    query()
