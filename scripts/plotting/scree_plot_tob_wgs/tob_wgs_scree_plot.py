"""
Make scree plot of TOB-WGS + NFE samples
(from HGDP + 1kg dataset) in order to identify
number of PCs to use in LM. Dependent on output from
`hail_batch/final_tob_batch_pca/generate_pca_no_outliers.py`
"""

import pandas as pd
import hail as hl
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.plotting import figure
from analysis_runner import bucket_path, output_path

EIGENVALUES = bucket_path('kat/pca/nfe/v0/eigenvalues.ht')


def query():  # pylint: disable=too-many-locals
    """Query script entry point."""

    eigenvalues = pd.read_csv(EIGENVALUES)

    # add in each PC component as a new column
    eigenvalues['PCs'] = list(eigenvalues.index + 1)
    # rename '0' column to 'variance'
    eigenvalues = eigenvalues.rename(columns={eigenvalues.columns[0]: 'variance'})
    # convert df to strings in order to plot
    eigenvalues = eigenvalues.applymap(str)

    # plot
    p = figure(
        x_range=eigenvalues.PCs,
        height=500,
        title='Variance explained',
        toolbar_location=None,
        tools='',
    )
    p.vbar(x=eigenvalues.PCs, top=eigenvalues.variance, width=0.9)
    p.xaxis.axis_label = 'PC'
    p.yaxis.axis_label = 'Percentage of variance explained'
    p.xgrid.grid_line_color = None
    p.y_range.start = 0
    # save plot
    plot_filename_html = output_path(f'scree_plot_tob_wgs.html', 'web')
    html = file_html(p, CDN, 'my plot')
    plot_filename_html = output_path(f'scree_plot.html', 'web')
    with hl.hadoop_open(plot_filename_html, 'w') as f:
        f.write(html)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
