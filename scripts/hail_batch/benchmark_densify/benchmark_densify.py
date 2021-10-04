"""Test cost and speed of densifying on 20 samples"""

import hail as hl
from analysis_runner import output_path

MT = 'gs://cpg-tob-wgs-test-tmp/joint-calling/v6-25/combiner/v6-25-raw.mt'
SPLIT_MT = 'gs://cpg-tob-wgs-test/mt/v5.1.mt/'


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(MT)
    split_mt = hl.read_matrix_table(SPLIT_MT)
    mt = mt.select_entries('LGT', 'END').select_cols().select_rows()
    split_mt = split_mt.select_entries('GT', 'END').select_cols().select_rows()
    mt_path = output_path(f'v6-25-raw-filtered.mt', 'tmp')
    mt.write(mt_path)
    split_mt_path = output_path(f'v5.1-filtered.mt', 'tmp')
    split_mt.write(split_mt_path)


if __name__ == '__main__':
    query()
