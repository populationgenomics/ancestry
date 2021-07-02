"""Create PCA plots for the combined TOB-WGS/SNP-chip data"""

from bokeh.io.export import get_screenshot_as_png
from bokeh.resources import CDN
from bokeh.embed import file_html
import pandas as pd
import hail as hl
import click
from analysis_runner import bucket_path, output_path

SCORES = bucket_path('tob_wgs_snp_chip_variant_pca/v6/scores.ht/')
EIGENVALUES = bucket_path('tob_wgs_snp_chip_variant_pca/v6/eigenvalues.ht')


@click.command()
def query():  # pylint: disable=too-many-locals
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    scores = hl.read_table(SCORES)
    scores = scores.annotate(
        cohort_sample_codes=hl.if_else(
            scores.s.contains('snp_chip'), 'snp_chip', 'tob_wgs'
        )
    )
    labels = scores.cohort_sample_codes
    hover_fields = dict([('s', scores.s)])

    # get percent variance explained
    eigenvalues = hl.import_table(EIGENVALUES)
    eigenvalues = eigenvalues.to_pandas()
    eigenvalues.columns = ['eigenvalue']
    eigenvalues = pd.to_numeric(eigenvalues.eigenvalue)
    variance = eigenvalues.divide(float(eigenvalues.sum())) * 100
    variance = variance.round(2)

    # Get number of PCs
    number_of_pcs = len(eigenvalues)

    for i in range(0, (number_of_pcs - 1)):
        pc1 = i
        pc2 = i + 1
        print(f'PC{pc1 + 1} vs PC{pc2 + 1}')
        p = hl.plot.scatter(
            scores.scores[pc1],
            scores.scores[pc2],
            label=labels,
            title='TOB-WGS + TOB SNP Chip',
            xlabel='PC' + str(pc1 + 1) + ' (' + str(variance[pc1]) + '%)',
            ylabel='PC' + str(pc2 + 1) + ' (' + str(variance[pc2]) + '%)',
            hover_fields=hover_fields,
        )
        plot_filename = output_path('pc' + str(pc2) + '.png', 'web')
        with hl.hadoop_open(plot_filename, 'wb') as f:
            get_screenshot_as_png(p).save(f, format='PNG')
        html = file_html(p, CDN, 'my plot')
        plot_filename_html = output_path('pc' + str(pc2) + '.html', 'web')
        with hl.hadoop_open(plot_filename_html, 'w') as f:
            f.write(html)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
