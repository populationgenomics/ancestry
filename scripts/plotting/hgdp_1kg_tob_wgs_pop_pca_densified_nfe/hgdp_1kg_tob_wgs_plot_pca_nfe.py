"""Create PCA plots for HGDP/1kG + TOB-WGS samples"""

import subprocess
from bokeh.models import CategoricalColorMapper
from bokeh.palettes import turbo  # pylint: disable=no-name-in-module
from bokeh.io.export import get_screenshot_as_png
from bokeh.plotting import output_file, save
import pandas as pd
import hail as hl
import click

HGDP1KG_TOBWGS = (
    'gs://cpg-tob-wgs-main/1kg_hgdp_densified_pca/v2/'
    'hgdp1kg_tobwgs_joined_all_samples.mt/'
)
SCORES = 'gs://cpg-tob-wgs-main/1kg_hgdp_densified_nfe/v0/scores.ht/'
EIGENVALUES = 'gs://cpg-tob-wgs-main/1kg_hgdp_densified_nfe/v0/eigenvalues.ht/'


@click.command()
@click.option('--output', help='GCS output path', required=True)
def query(output):  # pylint: disable=too-many-locals
    """Query script entry point."""

    hl.init()

    mt = hl.read_matrix_table(HGDP1KG_TOBWGS)
    # Get NFE samples only
    mt = mt.filter_cols(
        (mt.hgdp_1kg_metadata.population_inference.pop == 'nfe')
        | (mt.s.contains('TOB'))
    )
    scores = hl.read_table(SCORES)
    mt = mt.annotate_cols(scores=scores[mt.s].scores)
    mt = mt.annotate_cols(TOB_WGS=mt.s.contains('TOB'))

    # PCA plot must all come from the same object
    columns = mt.cols()
    pca_scores = columns.scores
    labels = columns.TOB_WGS
    hover_fields = dict([('s', columns.s)])

    # get percent variance explained
    eigenvalues = hl.import_table(EIGENVALUES)
    eigenvalues = eigenvalues.to_pandas()
    eigenvalues.columns = ['eigenvalue']
    eigenvalues = pd.to_numeric(eigenvalues.eigenvalue)
    variance = eigenvalues.divide(float(eigenvalues.sum())) * 100
    variance = variance.round(2)

    # Get number of PCs
    number_of_pcs = len(eigenvalues)

    print('Making PCA plots labelled by study')
    for i in range(0, (number_of_pcs - 1)):
        pc1 = i
        pc2 = i + 1
        print(f'PC{pc1 + 1} vs PC{pc2 + 1}')
        p = hl.plot.scatter(
            pca_scores[pc1],
            pca_scores[pc2],
            label=labels,
            title='TOB-WGS',
            xlabel='PC' + str(pc1 + 1) + ' (' + str(variance[pc1]) + '%)',
            ylabel='PC' + str(pc2 + 1) + ' (' + str(variance[pc2]) + '%)',
            hover_fields=hover_fields,
        )
        plot_filename = f'{output}/study_pc' + str(pc2) + '.png'
        with hl.hadoop_open(plot_filename, 'wb') as f:
            get_screenshot_as_png(p).save(f, format='PNG')
        plot_filename_html = 'study_pc' + str(pc2) + '.html'
        output_file(plot_filename_html)
        save(p)
        subprocess.run(['gsutil', 'cp', plot_filename_html, output], check=False)

    print('Making PCA plots labelled by the subpopulation')
    labels = columns.hgdp_1kg_metadata.labeled_subpop
    pops = list(set(labels.collect()))

    for i in range(0, (number_of_pcs - 1)):
        pc1 = i
        pc2 = i + 1
        print(f'PC{pc1 + 1} vs PC{pc2 + 1}')
        p = hl.plot.scatter(
            pca_scores[pc1],
            pca_scores[pc2],
            label=labels,
            title='Sub-Population',
            xlabel='PC' + str(pc1 + 1) + ' (' + str(variance[pc1]) + '%)',
            ylabel='PC' + str(pc2 + 1) + ' (' + str(variance[pc2]) + '%)',
            collect_all=True,
            colors=CategoricalColorMapper(palette=turbo(len(pops)), factors=pops),
        )
        plot_filename = f'{output}/subpopulation_pc' + str(pc2) + '.png'
        with hl.hadoop_open(plot_filename, 'wb') as f:
            get_screenshot_as_png(p).save(f, format='PNG')
        plot_filename_html = 'subpopulation_pc' + str(pc2) + '.html'
        output_file(plot_filename_html)
        save(p)
        subprocess.run(['gsutil', 'cp', plot_filename_html, output], check=False)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
