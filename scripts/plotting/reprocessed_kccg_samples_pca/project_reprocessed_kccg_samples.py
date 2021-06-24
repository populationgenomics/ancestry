"""
Perform PCA on samples reprocessed using
the KCCG pipeline.
"""

import subprocess
import click
import pandas as pd
import hail as hl
from hail.experimental import lgt_to_gt
from hail.experimental import pc_project
from bokeh.palettes import Dark2  # pylint: disable=no-name-in-module
from bokeh.io.export import get_screenshot_as_png
from bokeh.plotting import output_file, save, ColumnDataSource, figure
from bokeh.transform import factor_cmap


HGDP1KG_TOBWGS = (
    'gs://cpg-tob-wgs-main/1kg_hgdp_densified_pca/v2/'
    'hgdp1kg_tobwgs_joined_all_samples.mt/'
)
LOADINGS = 'gs://cpg-tob-wgs-main/1kg_hgdp_densified_nfe/v0/loadings.ht/'
REPROCESSED_1KG = 'gs://cpg-tob-wgs-test/pipeline_validation/kccg_gatk4/mt/v0.mt'
SCORES = 'gs://cpg-tob-wgs-main/1kg_hgdp_densified_nfe/v0/scores.ht/'
EIGENVALUES = 'gs://cpg-tob-wgs-main/1kg_hgdp_densified_nfe/v0/eigenvalues.ht/'


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
    reprocessed_samples = hl.read_matrix_table(REPROCESSED_1KG)
    reprocessed_samples = hl.experimental.densify(reprocessed_samples)
    reprocessed_samples = reprocessed_samples.annotate_entries(
        GT=lgt_to_gt(reprocessed_samples.LGT, reprocessed_samples.LA)
    )
    ht = pc_project(reprocessed_samples.GT, loadings.loadings, loadings.af)
    ht = ht.key_by(s=ht.s + '_reprocessed')
    pcs = hl.read_table(SCORES)
    union_scores = ht.union(pcs)
    union_scores = union_scores.annotate(
        original=(union_scores.s == 'HG01513')
        | (union_scores.s == 'HG02238')
        | (union_scores.s == 'NA12248')
        | (union_scores.s == 'NA20502')
        | (union_scores.s == 'NA20826'),
        reprocessed=union_scores.s.contains('reprocessed'),
    )
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
    eigenvalues = hl.import_table(EIGENVALUES)
    eigenvalues = eigenvalues.to_pandas()
    eigenvalues.columns = ['eigenvalue']
    eigenvalues = pd.to_numeric(eigenvalues.eigenvalue)
    variance = eigenvalues.divide(float(eigenvalues.sum())) * 100
    variance = variance.round(2)

    # plot
    labels = union_scores.cohort_sample_codes
    sample_names = union_scores.s
    cohort_sample_codes = list(set(labels.collect()))
    tooltips = [('labels', '@label'), ('samples', '@samples')]
    for i in range(0, 10):
        pc1 = i
        pc2 = i + 1
        plot_filename = (
            f'{output}/reprocessed_sample_projection_pc' + str(i + 1) + '.png'
        )
        if not hl.hadoop_exists(plot_filename):
            plot = figure(
                title='Reprocessed Sample Projection',
                x_axis_label='PC' + str(pc1 + 1) + ' (' + str(variance[pc1]) + '%)',
                y_axis_label='PC' + str(pc2 + 1) + ' (' + str(variance[pc1]) + '%)',
                tooltips=tooltips,
            )
            source = ColumnDataSource(
                dict(
                    x=union_scores.scores[pc1].collect(),
                    y=union_scores.scores[pc2].collect(),
                    label=labels.collect(),
                    samples=sample_names.collect(),
                )
            )
            plot.circle(
                'x',
                'y',
                alpha=0.5,
                source=source,
                size=8,
                color=factor_cmap(
                    'label', Dark2[len(cohort_sample_codes)], cohort_sample_codes
                ),
                legend_group='label',
            )
            plot.add_layout(plot.legend[0], 'left')
            with hl.hadoop_open(plot_filename, 'wb') as f:
                get_screenshot_as_png(plot).save(f, format='PNG')
            plot_filename_html = (
                'reprocessed_sample_projection_pc' + str(i + 1) + '.html'
            )
            output_file(plot_filename_html)
            save(plot)
            subprocess.run(['gsutil', 'cp', plot_filename_html, output], check=False)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
