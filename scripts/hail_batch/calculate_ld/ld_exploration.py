"""Calculate ld using the ld_matrix function"""

import hail as hl
from analysis_runner import bucket_path

LD_MATRIX = bucket_path('kat/v0/ld_matrix_full_filtered.ht/')
TOB_LOCUS_INFO = 'gs://cpg-tob-wgs-test/kat/v0/tob_locus_info.ht/'


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    ht = hl.read_table(LD_MATRIX)
    annotation = hl.read_table(TOB_LOCUS_INFO).key_by('row_idx')
    # add in global positions for i
    ht = ht.rename({'i': 'row_idx'}).key_by('row_idx')
    ht = ht.annotate(global_pos_i=annotation[ht.row_idx].global_position)
    # add in global positions for j
    ht = ht.rename({'row_idx': 'i', 'j': 'row_idx'}).key_by('row_idx')
    ht = ht.annotate(global_pos_j=annotation[ht.row_idx].global_position)
    # collapse ht by variant and save
    ht_agg = ht.key_by('global_pos_i').drop('i', 'row_idx').collect_by_key()
    ht_agg.write('gs://cpg-tob-wgs-test/kat/v0/ld_matrix_aggregated_table.ht')


if __name__ == '__main__':
    query()
