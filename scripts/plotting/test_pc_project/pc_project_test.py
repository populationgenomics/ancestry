"""
Perform PCA on three 1kG samples which have not been
reprocessed in order to test the pc project function.
"""

import click
import pandas as pd
import hail as hl
from hail.experimental import pc_project
from bokeh.models import CategoricalColorMapper
from bokeh.palettes import turbo  # pylint: disable=no-name-in-module
from bokeh.io.export import get_screenshot_as_png

HGDP1KG_TOBWGS = (
    'gs://cpg-tob-wgs-analysis/1kg_hgdp_tobwgs_pca/'
    'v1/hgdp1kg_tobwgs_joined_all_samples.mt'
)

LOADINGS = 'gs://cpg-tob-wgs-analysis/1kg_hgdp_nfe/v0/loadings.ht'
GNOMAD_HGDP_1KG_MT = (
    'gs://gcp-public-data--gnomad/release/3.1/mt/genomes/'
    'gnomad.genomes.v3.1.hgdp_1kg_subset_dense.mt'
)
SCORES = 'gs://cpg-tob-wgs-analysis/1kg_hgdp_nfe/v0/scores.ht'
EIGENVALUES = 'gs://cpg-tob-wgs-analysis/1kg_hgdp_nfe/v0/eigenvalues.csv'


@click.command()
@click.option('--output', help='GCS output path', required=True)
@click.option('--pop', help='Population to subset from the 1KG (e.g. afr, nfe)')
def query(output, pop):  # pylint: disable=too-many-locals
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(HGDP1KG_TOBWGS)
    if pop:
        # Get samples from the specified population only
        mt = mt.filter_cols(
            (mt.hgdp_1kg_metadata.population_inference.pop == pop.lower())
            | (mt.s.contains('TOB'))
        )
    else:
        mt = mt.filter_cols(mt.s.contains('TOB'))

    mt = mt.annotate_rows(af=hl.agg.mean(mt.GT.n_alt_alleles()) / 2)
    loadings = hl.read_table(LOADINGS)
    loadings = loadings.annotate(af=mt.rows()[loadings.key].af)
    hgdp_1kg_samples = hl.read_matrix_table(GNOMAD_HGDP_1KG_MT)
    samples_to_keep = hl.literal(['HG02238', 'NA12248', 'NA20502'])
    hgdp_1kg_samples = hgdp_1kg_samples.filter_cols(
        samples_to_keep.contains(hgdp_1kg_samples['s'])
    )
    ht = pc_project(hgdp_1kg_samples.GT, loadings.loadings, loadings.af)
    ht = ht.key_by(s=ht.s + '_reprocessed')
    pcs = hl.read_table(SCORES)
    union_scores = ht.union(pcs)
    union_scores = union_scores.annotate(
        original=samples_to_keep.contains(union_scores.s),
        reprocessed=union_scores.s.contains('reprocessed'),
    )
    # union_scores.annotate(original=hl.literal(names).contains(union_scores.s))
    expr = (
        hl.case()
        .when(
            (union_scores.original)
            & (
                union_scores.reprocessed  # pylint: disable=singleton-comparison
                == False  # noqa: E712
            ),
            'original',
        )
        .when(
            (union_scores.original == False)  # pylint: disable=singleton-comparison
            & (union_scores.reprocessed),
            'reprocessed',
        )
        .default('unedited')
    )
    union_scores = union_scores.annotate(cohort_sample_codes=expr)
    # get percentage of variance explained
    eigenvalues = pd.read_csv(EIGENVALUES)
    eigenvalues.columns = ['eigenvalue']
    variance = eigenvalues['eigenvalue'].divide(float(eigenvalues.sum())) * 100
    variance = variance.round(2)

    # plot
    labels = union_scores.cohort_sample_codes
    cohort_sample_codes = list(set(labels.collect()))
    for i in range(0, 10):
        pc1 = i
        pc2 = i + 1
        plot_filename = f'{output}/test_sample_projection_pc' + str(i + 1) + '.png'
        if not hl.hadoop_exists(plot_filename):
            plot = hl.plot.scatter(
                union_scores.scores[pc1],
                union_scores.scores[pc2],
                label=labels,
                title='Reprocessed Sample Projection',
                xlabel='PC' + str(pc1 + 1) + ' (' + str(variance[pc1]) + '%)',
                ylabel='PC' + str(pc2 + 1) + ' (' + str(variance[pc1]) + '%)',
                colors=CategoricalColorMapper(
                    palette=turbo(len(cohort_sample_codes)), factors=cohort_sample_codes
                ),
            )
            with hl.hadoop_open(plot_filename, 'wb') as f:
                get_screenshot_as_png(plot).save(f, format='PNG')


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
