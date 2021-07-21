"""
Project WGS data onto SNP-chip data
"""

import re
import hail as hl
from analysis_runner import bucket_path, output_path
from hail.experimental import pc_project
from hail.experimental import lgt_to_gt
from bokeh.plotting import ColumnDataSource, figure
from bokeh.palettes import Dark2  # pylint: disable=no-name-in-module
from bokeh.transform import factor_cmap
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.io.export import get_screenshot_as_png

SNP_CHIP = bucket_path(
    'tob_wgs_snp_chip_pca/increase_partitions/v0/snp_chip_100_partitions.mt'
)
# TOB_WGS = bucket_path('mt/v3-raw.mt')
TOB_WGS = bucket_path('mt/v4.mt')


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    snp_chip = hl.read_matrix_table(SNP_CHIP)
    tob_wgs = hl.read_matrix_table(TOB_WGS)
    tob_wgs = tob_wgs.head(1000000)
    tob_wgs = hl.experimental.densify(tob_wgs)
    tob_wgs = tob_wgs.annotate_entries(GT=lgt_to_gt(tob_wgs.LGT, tob_wgs.LA))

    snp_chip = snp_chip.semi_join_rows(tob_wgs.rows())
    # Perform PCA
    eigenvalues, scores, loadings = hl.hwe_normalized_pca(
        snp_chip.GT, compute_loadings=True, k=5
    )

    # make tob_wgs rows equivalent to the snp_chip rows
    tob_wgs = tob_wgs.semi_join_rows(snp_chip.rows())
    snp_chip = snp_chip.annotate_rows(af=hl.agg.mean(snp_chip.GT.n_alt_alleles()) / 2)
    loadings = loadings.annotate(af=snp_chip.rows()[loadings.key].af)
    # project WGS samples onto PCA
    ht = pc_project(tob_wgs.GT, loadings.loadings, loadings.af)
    scores = scores.key_by(s=scores.s + '_snp_chip')
    union_scores = ht.union(scores)
    variance = [(x / sum(eigenvalues) * 100) for x in eigenvalues]
    variance = [round(x, 2) for x in variance]

    # Get partner sample information
    sample_names = union_scores.s.collect()

    def sample_type(sample_name):
        if sample_name.endswith('snp_chip'):
            partner_name = re.sub('_snp_chip', '', sample_name)
            tech = 'snp'
        else:
            partner_name = sample_name + '_snp_chip'
            tech = 'wgs'

        if partner_name in sample_names:
            prefix = 'dual_'
        else:
            prefix = ''

        return prefix + tech

    # plot
    labels = list(map(sample_type, sample_names))
    cohort_sample_codes = list(set(labels))
    tooltips = [('labels', '@label'), ('samples', '@samples')]

    # Get number of PCs
    number_of_pcs = len(eigenvalues)

    for i in range(0, (number_of_pcs - 1)):
        pc1 = i
        pc2 = i + 1
        plot = figure(
            title='TOB-WGS + TOB SNP Chip',
            x_axis_label=f'PC{pc1 + 1} ({variance[pc1]})%)',
            y_axis_label=f'PC{pc2 + 1} ({variance[pc2]}%)',
            tooltips=tooltips,
        )
        source = ColumnDataSource(
            dict(
                x=union_scores.scores[pc1].collect(),
                y=union_scores.scores[pc2].collect(),
                label=labels,
                samples=sample_names,
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
        plot_filename = output_path(f'pc{pc2}.png', 'web')
        with hl.hadoop_open(plot_filename, 'wb') as f:
            get_screenshot_as_png(plot).save(f, format='PNG')
        html = file_html(plot, CDN, 'my plot')
        plot_filename_html = output_path(f'pc{pc2}.html', 'web')
        with hl.hadoop_open(plot_filename_html, 'w') as f:
            f.write(html)


if __name__ == '__main__':
    query()
