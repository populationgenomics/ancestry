"""Calculate ld using the ld_matrix function"""

import hail as hl
from bokeh.resources import CDN
from bokeh.embed import file_html
from analysis_runner import output_path

# from analysis_runner import bucket_path

# LD_MATRIX = bucket_path('kat/v0/ld_matrix_full_filtered02.ht/')
# TOB_LOCUS_INFO = 'gs://cpg-tob-wgs-test/kat/v0/tob_locus_info.ht/'


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    # ht = hl.read_table(LD_MATRIX)
    # annotation = hl.read_table(TOB_LOCUS_INFO).key_by('row_idx')
    # # add in global positions for i
    # ht = ht.rename({'i': 'row_idx'}).key_by('row_idx')
    # ht = ht.annotate(global_pos_i=annotation[ht.row_idx].global_position)
    # # add in global positions for j
    # ht = ht.rename({'row_idx': 'i', 'j': 'row_idx'}).key_by('row_idx')
    # ht = ht.annotate(global_pos_j=annotation[ht.row_idx].global_position)
    # # collapse ht by variant and save
    # ht_agg = ht.key_by('global_pos_i').drop('i', 'row_idx').collect_by_key()
    # ht_agg.write('gs://cpg-tob-wgs-test/kat/v0/ld_matrix_aggregated_table02.ht')
    ht_agg = hl.read_table(
        'gs://cpg-tob-wgs-test/kat/v0/ld_matrix_aggregated_table02.ht'
    )
    # min_val = ht_agg.aggregate(hl.agg.min(hl.len(ht_agg.values)))
    # max_val = ht_agg.aggregate(hl.agg.max(hl.len(ht_agg.values)))
    # median = ht_agg.aggregate(hl.agg.approx_quantiles(hl.len(ht_agg.values), 0.5))
    # print(f'The min is {min_val}, the max is {max_val}, and the median is {median}')
    ht_agg = ht_agg.annotate(n_ld=hl.len(ht_agg.values))
    plot = hl.plot.scatter(
        ht_agg.global_pos_i,
        ht_agg.n_ld,
        xlabel='Global Position',
        ylabel='ld variants (n)',
        title='Linked variants across the genome',
    )
    plot_filename_html = output_path(f'ld_variants_across_genome.html', 'web')
    html = file_html(plot, CDN, 'my plot')
    plot_filename_html = output_path(f'scree_plot.html', 'web')
    with hl.hadoop_open(plot_filename_html, 'w') as f:
        f.write(html)


if __name__ == '__main__':
    query()
