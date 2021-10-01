"""Test cost and speed of densifying on 20 samples"""

import hail as hl
from analysis_runner import output_path

SPARSE_MT = 'gs://cpg-tob-wgs-test-tmp/densify_benchmark/v5/sparse_tob_wgs.mt'
DENSIFIED_MT = 'gs://cpg-tob-wgs-test-tmp/densify_benchmark/v5/densified_tob_wgs.mt'


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    sparse = hl.read_matrix_table(SPARSE_MT)
    dense = hl.read_matrix_table(DENSIFIED_MT)
    # save entries field only
    dense_entries = dense.entries()
    dense_entries_path = output_path(f'dense_entries.ht', 'tmp')
    dense_entries.write(dense_entries_path)
    sparse_entries = sparse.entries()
    sparse_entries_path = output_path(f'sparse_entries.ht', 'tmp')
    sparse_entries.write(sparse_entries_path)
    # remove entries field and save
    dense = dense.select_entries()
    sparse = sparse.select_entries()
    dense_path = output_path(f'dense_tob_wgs_no_entries.mt', 'tmp')
    sparse_path = output_path(f'sparse_tob_wgs_no_entries.mt', 'tmp')
    dense.write(dense_path)
    sparse.write(sparse_path)


if __name__ == '__main__':
    query()
