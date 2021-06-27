"""
Project snp-chip samples onto
combined hgdp/1kg + tob-wgs dataset.
"""

import subprocess
import click
import pandas as pd
import hail as hl
from hail.experimental import pc_project
from bokeh.palettes import Dark2  # pylint: disable=no-name-in-module
from bokeh.io.export import get_screenshot_as_png
from bokeh.plotting import output_file, save, ColumnDataSource, figure
from bokeh.transform import factor_cmap


SNP_CHIP = 'gs://cpg-tob-wgs-test/snpchip/v1/snpchip_grch38.mt/'
# HGDP1KG_TOBWGS = (
#     'gs://cpg-tob-wgs-main/1kg_hgdp_densified_pca/v2/'
#     'hgdp1kg_tobwgs_joined_all_samples.mt/'
# )
# LOADINGS = 'gs://cpg-tob-wgs-main/1kg_hgdp_densified_nfe/v0/loadings.ht/'
# SCORES = 'gs://cpg-tob-wgs-main/1kg_hgdp_densified_nfe/v0/scores.ht/'
# EIGENVALUES = 'gs://cpg-tob-wgs-main/1kg_hgdp_densified_nfe/v0/eigenvalues.ht/'
HGDP1KG_TOBWGS = (
    'gs://cpg-tob-wgs-test/1kg_hgdp_tobwgs_pca/v1/hgdp1kg_tobwgs_joined_all_samples.mt/'
)
SCORES = 'gs://cpg-tob-wgs-test/1kg_hgdp_densify/v15/scores.ht/'
EIGENVALUES = 'gs://cpg-tob-wgs-test/1kg_hgdp_tobwgs_pca/v1/eigenvalues.ht'
LOADINGS = 'gs://cpg-tob-wgs-test/1kg_hgdp_densify/v15/loadings.ht/'


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
    tob_wgs_snp_chip = hl.read_matrix_table(SNP_CHIP).key_rows_by('locus', 'alleles')
    ht = pc_project(tob_wgs_snp_chip.GT, loadings.loadings, loadings.af)
    ht = ht.key_by(s=ht.s + '_SNP_CHIP')
    pcs = hl.read_table(SCORES)
    union_scores = ht.union(pcs)
    union_scores = union_scores.annotate(
        snp_chip=(union_scores.s.contains('_SNP_CHIP')),
        tob_wgs=(union_scores.s.contains('_SNP_CHIP') | union_scores.s.contains('TOB')),
    )
    expr = (
        hl.case()
        .when(
            (union_scores.snp_chip),
            'snp_chip',
        )
        .when(
            (
                union_scores.snp_chip  # noqa: E501; pylint: disable=singleton-comparison;
                == False  # noqa: E712
            )
            & (union_scores.tob_wgs),
            'tob_wgs',
        )
        .default('hgdp_1kg')
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
    number_of_pcs = len(eigenvalues)
    union_scores = union_scores.persist()
    for i in range(0, (number_of_pcs - 1)):
        pc1 = i
        pc2 = i + 1
        plot_filename = (
            f'{output}/reprocessed_sample_projection_pc' + str(i + 1) + '.png'
        )
        if not hl.hadoop_exists(plot_filename):
            plot = figure(
                title='SNP-Chip Sample Projection',
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
            plot_filename_html = 'snp_chip_sample_projection_pc' + str(i + 1) + '.html'
            output_file(plot_filename_html)
            save(plot)
            subprocess.run(['gsutil', 'cp', plot_filename_html, output], check=False)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
