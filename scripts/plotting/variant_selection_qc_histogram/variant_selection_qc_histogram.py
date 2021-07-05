"""QC of newly-selected variants"""

import click
import hail as hl
import numpy as np
import pandas as pd
from analysis_runner import output_path

# from analysis_runner import bucket_path, output_path
from bokeh.plotting import figure
from bokeh.io.export import get_screenshot_as_png
from bokeh.resources import CDN
from bokeh.embed import file_html

# FILTERED_VARIANTS = bucket_path(
#     'tob_wgs_hgdp_1kg_variant_selection/v8/'
#     'tob_wgs_hgdp_1kg_filtered_variants.mt'
# )

FILTERED_VARIANTS = (
    'gs://cpg-tob-wgs-test/1kg_hgdp_tobwgs_pca/'
    'v1/hgdp1kg_tobwgs_joined_all_samples.mt/'
)


@click.command()
def query():  # pylint: disable=too-many-locals
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(FILTERED_VARIANTS)
    mt = mt.head(100)
    mt = hl.variant_qc(mt)
    print(mt.count_rows)

    # Plot the allele frequency
    fig = figure(
        title='Variant AF',
        x_axis_label='Allele Frequency',
        y_axis_label='Frequency (%)',
    )
    variant_af = mt.variant_qc.AF[1].collect()
    af_count, edges = np.histogram(
        variant_af, weights=np.ones(len(variant_af)) / len(variant_af)
    )
    variant_af_count = pd.DataFrame(
        {'variant_af_count': af_count, 'left': edges[:-1], 'right': edges[1:]}
    )
    fig.quad(
        bottom=0,
        top=variant_af_count['variant_af_count'],
        left=variant_af_count['left'],
        right=variant_af_count['right'],
        fill_color='blue',
        line_color='black',
    )
    # Add in the cumulative distribution
    cumulative_af = np.cumsum(af_count)
    fig.line(
        x=variant_af_count['right'],
        y=cumulative_af,
        color='gray',
        line_width=1,
        legend='Cum dist',
    )
    fig.legend.location = 'top_left'
    fig_filename = output_path('variant_selection_histogram.png', 'web')
    with hl.hadoop_open(fig_filename, 'wb') as f:
        get_screenshot_as_png(fig).save(f, format='PNG')
    html = file_html(fig, CDN, 'my plot')
    fig_filename_html = output_path('variant_selection_histogram.html', 'web')
    with hl.hadoop_open(fig_filename_html, 'w') as f:
        f.write(html)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
