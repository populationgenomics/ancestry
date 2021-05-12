"""Create PCA plots for HGDP/1kG + TOB-WGS samples"""

from bokeh.models import CategoricalColorMapper
from bokeh.palettes import turbo  # pylint: disable=no-name-in-module
import pandas as pd
import hail as hl
from hail.plot import show
import click


HGDP1KG_TOBWGS = (
    'gs://cpg-tob-wgs-analysis/1kg_hgdp_tobwgs_pca/v1/'
    'hgdp1kg_tobwgs_joined_all_samples.mt'
)
SCORES = 'gs://cpg-tob-wgs-analysis/1kg_hgdp_tobwgs_pca/v1/scores.ht/'
EIGENVALUES = 'gs://cpg-tob-wgs-analysis/1kg_hgdp_tobwgs_pca/v1/eigenvalues.csv'


@click.command()
@click.option('--number-of-pcs', 'number_of_pcs', help='Number of PCS', default=19)
def main(number_of_pcs: int):  # pylint: disable=too-many-locals
    hl.init()
    
    mt = hl.read_matrix_table(HGDP1KG_TOBWGS)
    scores = hl.read_table(SCORES)
    mt = mt.annotate_cols(scores=scores[mt.s].scores)
    mt = mt.annotate_cols(TOB_WGS=mt.s.contains('TOB'))
    
    # PCA plot must all come from the same object
    columns = mt.cols()
    pca_scores = columns.scores
    labels = columns.TOB_WGS
    
    # get percent variance explained
    eigenvalues = pd.read_csv(EIGENVALUES)
    eigenvalues.columns = ['eigenvalue']
    variance = eigenvalues['eigenvalue'].divide(float(eigenvalues.sum())) * 100
    variance = variance.round(2)
    
    print('Making PCA plots labelled by the study ID')
    for i in range(0, number_of_pcs):
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
        )
        show(p)
    
    print('Making PCA plots labelled by the continental population')
    labels = columns.hgdp_1kg_metadata.population_inference.pop
    pops = list(set(labels.collect()))
    hover_fields = dict([('s', columns.s)])
    
    for i in range(0, number_of_pcs):
        pc1 = i
        pc2 = i + 1
        print(f'PC{pc1 + 1} vs PC{pc2 + 1}')
        p = hl.plot.scatter(
            pca_scores[pc1],
            pca_scores[pc2],
            label=labels,
            title='Continental Population',
            xlabel='PC' + str(pc1 + 1) + ' (' + str(variance[pc1]) + '%)',
            ylabel='PC' + str(pc2 + 1) + ' (' + str(variance[pc2]) + '%)',
            collect_all=True,
            colors=CategoricalColorMapper(palette=turbo(len(pops)), factors=pops),
            hover_fields=hover_fields,
        )
        show(p)
    
    print('Making PCA plots labelled by the subpopulation')
    labels = columns.hgdp_1kg_metadata.labeled_subpop
    pops = list(set(labels.collect()))
    
    for i in range(0, number_of_pcs):
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
        show(p)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
