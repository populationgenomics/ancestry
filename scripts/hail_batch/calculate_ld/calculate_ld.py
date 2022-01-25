"""Calculate ld using the ld_matrix function"""

import hail as hl
from analysis_runner import bucket_path

TOB_WGS = bucket_path('mt/v7.mt/')


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    tob_wgs = hl.read_matrix_table(TOB_WGS)
    # add in global position and stats for each locus
    tob_wgs = tob_wgs.annotate_rows(
        global_position=hl.locus(
            tob_wgs.locus.contig, tob_wgs.locus.position
        ).global_position(),
        stats=hl.agg.stats(tob_wgs.GT.n_alt_alleles()),
    )
    tob_wgs = tob_wgs.filter_rows(tob_wgs.stats.stdev != 0)
    # add row index to be able to remap
    tob_wgs = tob_wgs.add_row_index()
    # turn tob matrix into table and save
    tob_wgs = tob_wgs.key_rows_by('locus', 'alleles', 'row_idx', 'global_position')
    tob_wgs = tob_wgs.select_rows().select_globals()
    tob_locus_info = tob_wgs.rows()
    tob_locus_info.write('gs://cpg-tob-wgs-test/kat/v0/tob_locus_info.ht')
    ld = hl.ld_matrix(tob_wgs.GT.n_alt_alleles(), tob_wgs.locus, radius=2e6)
    table = ld.entries()
    # filter out entries with an LD score less than 0.1
    table = table.filter(table.entry > 0.1)
    # save table
    table.write('gs://cpg-tob-wgs-test/kat/v0/ld_matrix_full_filtered.ht')


if __name__ == '__main__':
    query()
