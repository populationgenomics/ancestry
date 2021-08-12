"""
Create PCA plots for nfe samples
from the HGDP/1kG + TOB-WGS datasets,
removing outliers.
"""

import hail as hl
import pandas as pd
from analysis_runner import bucket_path, output_path
from bokeh.io.export import get_screenshot_as_png
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.plotting import ColumnDataSource, figure
from bokeh.transform import factor_cmap
from bokeh.palettes import turbo  # pylint: disable=no-name-in-module

HGDP1KG_TOBWGS = bucket_path(
    '1kg_hgdp_densified_pca_new_variants/v0/hgdp1kg_tobwgs_joined_all_samples.mt'
)
SCORES = bucket_path('tob_wgs_hgdp_1kg_nfe_pca_new_variants/v4/scores.ht/')
EIGENVALUES = bucket_path('tob_wgs_hgdp_1kg_nfe_pca_new_variants/v4/eigenvalues.ht')


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(HGDP1KG_TOBWGS)
    scores = hl.read_table(SCORES)
    # Get filtered NFE samples and plot whether they are included in gnomAD
    mt = mt.filter_cols(
        (mt.hgdp_1kg_metadata.population_inference.pop == 'nfe')
        | (mt.s.contains('TOB'))
    )
    filtered_samples = mt.filter_cols(hl.is_defined(scores[mt.col_key]), keep=False)
    filt_samples_in_gnomad = filtered_samples.hgdp_1kg_metadata.gnomad_release

    # save as html
    html = pd.DataFrame(
        {
            'sample_name': filtered_samples.s.collect(),
            'in_gnomad': filt_samples_in_gnomad.collect(),
        }
    ).to_html()
    plot_filename_html = output_path(f'filtered_samples_in_gnomad.html', 'web')
    with hl.hadoop_open(plot_filename_html, 'w') as f:
        f.write(html)

    # Filter outliers and related samples
    mt = mt.semi_join_cols(scores)
    mt = mt.annotate_cols(scores=scores[mt.s].scores)
    mt = mt.annotate_cols(study=hl.if_else(mt.s.contains('TOB'), 'TOB-WGS', 'HGDP-1kG'))

    # PCA plot must all come from the same object
    columns = mt.cols()
    pca_scores = columns.scores
    labels = columns.study
    sample_names = columns.s
    cohort_sample_codes = list(set(labels.collect()))
    tooltips = [('labels', '@label'), ('samples', '@samples')]

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
        plot = figure(
            title='TOB-WGS + HGDP/1kG Dataset',
            x_axis_label=f'PC{pc1 + 1} ({variance[pc1]})%)',
            y_axis_label=f'PC{pc2 + 1} ({variance[pc2]})%)',
            tooltips=tooltips,
        )
        source = ColumnDataSource(
            dict(
                x=pca_scores[pc1].collect(),
                y=pca_scores[pc2].collect(),
                label=labels.collect(),
                samples=sample_names.collect(),
            )
        )
        plot.circle(
            'x',
            'y',
            alpha=0.5,
            source=source,
            size=4,
            color=factor_cmap('label', ['#1b9e77', '#d95f02'], cohort_sample_codes),
            legend_group='label',
        )
        plot.add_layout(plot.legend[0], 'left')
        plot_filename = output_path(f'study_pc{pc2}.png', 'web')
        with hl.hadoop_open(plot_filename, 'wb') as f:
            get_screenshot_as_png(plot).save(f, format='PNG')
        html = file_html(plot, CDN, 'my plot')
        plot_filename_html = output_path(f'study_pc{pc2}.html', 'web')
        with hl.hadoop_open(plot_filename_html, 'w') as f:
            f.write(html)

    print('Making PCA plots labelled by the subpopulation')
    labels = columns.hgdp_1kg_metadata.labeled_subpop.collect()
    labels = ['TOB-WGS' if x is None else x for x in labels]
    subpopulation = list(set(labels))
    # change ordering of subpopulations
    # so TOB-WGS is at the end and glyphs appear on top
    subpopulation.append(subpopulation.pop(subpopulation.index('TOB-WGS')))
    tooltips = [('labels', '@label'), ('samples', '@samples')]

    for i in range(0, (number_of_pcs - 1)):
        pc1 = i
        pc2 = i + 1
        print(f'PC{pc1 + 1} vs PC{pc2 + 1}')
        plot = figure(
            title='Subpopulation',
            x_axis_label=f'PC{pc1 + 1} ({variance[pc1]}%)',
            y_axis_label=f'PC{pc2 + 1} ({variance[pc2]}%)',
            tooltips=tooltips,
        )
        source = ColumnDataSource(
            dict(
                x=pca_scores[pc1].collect(),
                y=pca_scores[pc2].collect(),
                label=labels,
                samples=sample_names.collect(),
            )
        )
        plot.circle(
            'x',
            'y',
            alpha=0.5,
            source=source,
            size=4,
            color=factor_cmap('label', turbo(len(subpopulation)), subpopulation),
            legend_group='label',
        )
        plot.add_layout(plot.legend[0], 'left')
        plot_filename = output_path(f'subpopulation_pc{pc2}.png', 'web')
        with hl.hadoop_open(plot_filename, 'wb') as f:
            get_screenshot_as_png(plot).save(f, format='PNG')
        html = file_html(plot, CDN, 'my plot')
        plot_filename_html = output_path(f'subpopulation_pc{pc2}.html', 'web')
        with hl.hadoop_open(plot_filename_html, 'w') as f:
            f.write(html)


if __name__ == '__main__':
    query()
