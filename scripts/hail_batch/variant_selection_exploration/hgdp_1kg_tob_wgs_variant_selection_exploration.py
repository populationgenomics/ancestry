"""Pipeline for choosing new variants for HGDP/1kG + TOB-WGS data"""

import click
import hail as hl
from hail.experimental import lgt_to_gt
from bokeh.io.export import get_screenshot_as_png

GNOMAD_HGDP_1KG_MT = (
    'gs://gcp-public-data--gnomad/release/3.1/mt/genomes/'
    'gnomad.genomes.v3.1.hgdp_1kg_subset_dense.mt'
)

TOB_WGS = 'gs://cpg-tob-wgs-main/mt/v2-raw.mt/'


@click.command()
@click.option('--output', help='GCS output path', required=True)
def query(output):  # pylint: disable=too-many-locals
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    hgdp_1kg = hl.read_matrix_table(GNOMAD_HGDP_1KG_MT)
    tob_wgs = hl.read_matrix_table(TOB_WGS).key_rows_by('locus', 'alleles')

    # filter to loci that are contained in both matrix tables after densifying
    tob_wgs = hl.experimental.densify(tob_wgs)
    hgdp_1kg = hgdp_1kg.filter_rows(
        hl.is_defined(tob_wgs.index_rows(hgdp_1kg['locus'], hgdp_1kg['alleles']))
    )
    tob_wgs = tob_wgs.semi_join_rows(hgdp_1kg.rows())

    # Entries and columns must be identical
    tob_wgs_select = tob_wgs.select_entries(GT=lgt_to_gt(tob_wgs.LGT, tob_wgs.LA))
    hgdp_1kg_select = hgdp_1kg.select_entries(hgdp_1kg.GT)
    hgdp_1kg_select = hgdp_1kg_select.select_cols()
    # Join datasets
    hgdp1kg_tobwgs_joined = hgdp_1kg_select.union_cols(tob_wgs_select)
    # Add in metadata information
    hgdp_1kg_metadata = hgdp_1kg.cols()
    hgdp1kg_tobwgs_joined = hgdp1kg_tobwgs_joined.annotate_cols(
        hgdp_1kg_metadata=hgdp_1kg_metadata[hgdp1kg_tobwgs_joined.s]
    )

    # choose variants based off of gnomAD v3 parameters
    hgdp1kg_tobwgs_joined = hl.variant_qc(hgdp1kg_tobwgs_joined)
    hgdp1kg_tobwgs_joined = hgdp1kg_tobwgs_joined.annotate_rows(
        IB=hl.agg.inbreeding(
            hgdp1kg_tobwgs_joined.GT, hgdp1kg_tobwgs_joined.variant_qc.AF[1]
        )
    )
    hgdp1kg_tobwgs_joined = hgdp1kg_tobwgs_joined.filter_rows(
        (hl.len(hgdp1kg_tobwgs_joined.alleles) == 2)
        & (hgdp1kg_tobwgs_joined.locus.in_autosome())
        & (hgdp1kg_tobwgs_joined.variant_qc.AF[1] > 0.001)
        & (hgdp1kg_tobwgs_joined.variant_qc.call_rate > 0.99)
        & (hgdp1kg_tobwgs_joined.IB.f_stat > -0.25)
    )

    histogram_plot = hl.plot.histogram(
        hgdp1kg_tobwgs_joined.variant_qc.AF[1],
        legend='Allele Frequency',
    )
    plot_filename = f'{output}/histogram_plot.png'
    with hl.hadoop_open(plot_filename, 'wb') as f:
        get_screenshot_as_png(histogram_plot).save(f, format='PNG')


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
