"""Plot PCA loadings for HGDP/1kG + TOB-WGS samples"""

from bokeh.models import CategoricalColorMapper, HoverTool
from bokeh.io.export import get_screenshot_as_png
from bokeh.plotting import figure
from bokeh.embed import file_html
from bokeh.resources import CDN
from analysis_runner import bucket_path, output_path
import hail as hl

LOADINGS = bucket_path('tob_wgs_hgdp_1kg_nfe_pca_new_variants/v4/scores.ht/')


def manhattan_loadings(
    pvals,
    locus=None,
    title=None,
    size=4,
    hover_fields=None,
    collect_all=False,
    n_divisions=500,
):

    """modify hail manhattan plot"""
    palette = [
        '#1f77b4',
        '#ff7f0e',
        '#2ca02c',
        '#d62728',
        '#9467bd',
        '#8c564b',
        '#e377c2',
        '#7f7f7f',
        '#bcbd22',
        '#17becf',
    ]
    if locus is None:
        locus = pvals._indices.source.locus  # pylint: disable=protected-access

    ref = locus.dtype.reference_genome

    if hover_fields is None:
        hover_fields = {}

    hover_fields['locus'] = hl.str(locus)

    source_pd = (
        hl.plot.plots._collect_scatter_plot_data(  # pylint: disable=protected-access
            ('_global_locus', locus.global_position()),
            ('_pval', pvals),
            fields=hover_fields,
            n_divisions=None if collect_all else n_divisions,
        )
    )
    source_pd['p_value'] = source_pd['_pval']
    source_pd['_contig'] = [locus.split(':')[0] for locus in source_pd['locus']]

    observed_contigs = set(source_pd['_contig'])
    observed_contigs = [
        contig for contig in ref.contigs.copy() if contig in observed_contigs
    ]

    contig_ticks = [
        ref._contig_global_position(contig)  # pylint: disable=protected-access
        + ref.contig_length(contig) // 2
        for contig in observed_contigs
    ]
    color_mapper = CategoricalColorMapper(
        factors=ref.contigs, palette=palette[:2] * int((len(ref.contigs) + 1) / 2)
    )

    p = figure(
        title=title, x_axis_label='Chromosome', y_axis_label='Loadings', width=1000
    )
    (
        p,
        _,
        legend,
        _,
        _,
        _,
    ) = hl.plot.plots._get_scatter_plot_elements(  # pylint: disable=protected-access
        p,
        source_pd,
        x_col='_global_locus',
        y_col='_pval',
        label_cols=['_contig'],
        colors={'_contig': color_mapper},
        size=size,
    )
    legend.visible = False
    p.xaxis.ticker = contig_ticks
    p.xaxis.major_label_overrides = dict(zip(contig_ticks, observed_contigs))
    p.select_one(HoverTool).tooltips = [
        t for t in p.select_one(HoverTool).tooltips if not t[0].startswith('_')
    ]

    return p


def query():  # pylint: disable=too-many-locals
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    loadings_ht = hl.read_table(LOADINGS)
    number_of_pcs = hl.len(loadings_ht.loadings).take(1)[0]
    for i in range(0, (number_of_pcs)):
        pc = i + 1
        p = manhattan_loadings(
            pvals=hl.abs(loadings_ht.loadings[i]),
            locus=loadings_ht.locus,
            title=f'Loadings of PC{pc}',
            collect_all=True,
        )
        plot_filename = output_path(f'loadings_manhattan_plot_pc{pc}.png', 'web')
        with hl.hadoop_open(plot_filename, 'wb') as f:
            get_screenshot_as_png(p).save(f, format='PNG')
        html = file_html(p, CDN, 'my plot')
        plot_filename_html = output_path(f'loadings_pc{pc}.html', 'web')
        with hl.hadoop_open(plot_filename_html, 'w') as f:
            f.write(html)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
