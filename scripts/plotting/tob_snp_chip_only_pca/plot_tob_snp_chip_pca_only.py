"""Create PCA plots for the combined TOB-WGS/SNP-chip data"""

from bokeh.io.export import get_screenshot_as_png
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.transform import factor_cmap
from bokeh.plotting import ColumnDataSource, figure
import pandas as pd
import hail as hl
import click
from analysis_runner import bucket_path, output_path

SCORES = bucket_path('tob_snp_chip_pca/v0/scores.ht')
EIGENVALUES = bucket_path('tob_snp_chip_pca/v0/eigenvalues.ht')
TOB_WGS = bucket_path('mt/v3-raw.mt')


@click.command()
def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    scores = hl.read_table(SCORES)
    tob_wgs = hl.read_matrix_table(TOB_WGS)
    snp_chip_names = scores.s.collect()
    wgs_names = tob_wgs.s.collect()

    def sample_type(sample_name):
        return 'dual_sample' if sample_name in wgs_names else 'snp_chip_only'

    labels = list(map(sample_type, snp_chip_names))

    # get percent variance explained
    eigenvalues = hl.import_table(EIGENVALUES)
    eigenvalues = eigenvalues.to_pandas()
    eigenvalues.columns = ['eigenvalue']
    eigenvalues = pd.to_numeric(eigenvalues.eigenvalue)
    variance = eigenvalues.divide(float(eigenvalues.sum())) * 100
    variance = variance.round(2)

    # Get number of PCs
    number_of_pcs = len(eigenvalues)

    # plot
    cohort_sample_codes = list(set(labels))
    tooltips = [('labels', '@label'), ('samples', '@samples')]
    for i in range(0, (number_of_pcs - 1)):
        pc1 = i
        pc2 = i + 1
        plot = figure(
            title='Reprocessed Sample Projection',
            x_axis_label=f'PC{pc1 + 1} ({variance[pc1]})%)',
            y_axis_label=f'PC{pc2 + 1} ({variance[pc2]}%)',
            tooltips=tooltips,
        )
        source = ColumnDataSource(
            dict(
                x=scores.scores[pc1].collect(),
                y=scores.scores[pc2].collect(),
                label=labels,
                samples=snp_chip_names,
            )
        )
        plot.circle(
            'x',
            'y',
            alpha=0.5,
            source=source,
            size=8,
            color=factor_cmap('label', ['#1b9e77', '#d95f02'], cohort_sample_codes),
            legend_group='label',
        )
        plot.add_layout(plot.legend[0], 'left')
        plot_filename = output_path(f'pc{pc2}.png', 'web')
        with hl.hadoop_open(plot_filename, 'wb') as f:
            get_screenshot_as_png(plot).save(f, format='PNG')
        html = file_html(plot, CDN, 'my plot')
        plot_filename_html = output_path(f'pc{pc2}.html', 'web')
        with hl.hadoop_open(plot_filename_html, 'w') as f:
            f.write(html)


if __name__ == '__main__':
    query()
