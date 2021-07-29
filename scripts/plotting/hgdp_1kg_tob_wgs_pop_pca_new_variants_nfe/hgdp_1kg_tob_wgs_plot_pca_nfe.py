"""
Create PCA plots for nfe samples
from the HGDP/1kG + TOB-WGS datasets
"""

import hail as hl
import pandas as pd
from analysis_runner import bucket_path, output_path
from bokeh.io.export import get_screenshot_as_png
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.palettes import turbo  # pylint: disable=no-name-in-module
from bokeh.models import CategoricalColorMapper

HGDP1KG_TOBWGS = bucket_path(
    '1kg_hgdp_densified_pca_new_variants/v0/hgdp1kg_tobwgs_joined_all_samples.mt'
)
SCORES = bucket_path('1kg_hgdp_densified_nfe_new_variants/v0/scores.ht')
EIGENVALUES = bucket_path('1kg_hgdp_densified_nfe_new_variants/v0/eigenvalues.ht')


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

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
            collect_all=True,
            hover_fields=hover_fields,
        )
        plot_filename = output_path(f'study_pc{pc2}.png', 'web')
        with hl.hadoop_open(plot_filename, 'wb') as f:
            get_screenshot_as_png(p).save(f, format='PNG')
        html = file_html(p, CDN, 'my plot')
        plot_filename_html = output_path(f'study_pc{pc2}.html', 'web')
        with hl.hadoop_open(plot_filename_html, 'w') as f:
            f.write(html)

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
            title='Subpopulation',
            xlabel='PC' + str(pc1 + 1) + ' (' + str(variance[pc1]) + '%)',
            ylabel='PC' + str(pc2 + 1) + ' (' + str(variance[pc2]) + '%)',
            collect_all=True,
            colors=CategoricalColorMapper(palette=turbo(len(pops)), factors=pops),
        )
        plot_filename = output_path(f'subpopulation_pc{pc2}.png', 'web')
        with hl.hadoop_open(plot_filename, 'wb') as f:
            get_screenshot_as_png(p).save(f, format='PNG')
        html = file_html(p, CDN, 'my plot')
        plot_filename_html = output_path(f'subpopulation_pc{pc2}.html', 'web')
        with hl.hadoop_open(plot_filename_html, 'w') as f:
            f.write(html)


if __name__ == '__main__':
    query()
