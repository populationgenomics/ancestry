"""
Densify TOB-WGS data.
"""

import click
import hail as hl
from analysis_runner import bucket_path, output_path


NEW_VARIANTS = bucket_path(
    'tob_wgs_hgdp_1kg_variant_selection/v8/tob_wgs_hgdp_1kg_filtered_variants.mt/'
)
GNOMAD_HGDP_1KG_MT = (
    'gs://gcp-public-data--gnomad/release/3.1/mt/genomes/'
    'gnomad.genomes.v3.1.hgdp_1kg_subset_dense.mt'
)
TOB_WGS = bucket_path('mt/v4.mt//')


@click.command()
def query():  # pylint: disable=too-many-locals
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    hgdp_1kg = hl.read_matrix_table(GNOMAD_HGDP_1KG_MT)
    tob_wgs = hl.read_matrix_table(TOB_WGS).key_rows_by('locus', 'alleles')
    new_variants = hl.read_matrix_table(NEW_VARIANTS).key_rows_by('locus', 'alleles')

    # filter to loci that are contained in both tables and the loadings after densifying
    tob_wgs = hl.experimental.densify(tob_wgs)
    hgdp_1kg = hgdp_1kg.filter_rows(
        hl.is_defined(new_variants.index_rows(hgdp_1kg['locus'], hgdp_1kg['alleles']))
        & hl.is_defined(tob_wgs.index_rows(hgdp_1kg['locus'], hgdp_1kg['alleles']))
    )
    tob_wgs = tob_wgs.semi_join_rows(hgdp_1kg.rows())
    tob_wgs = tob_wgs.cache()
    print(tob_wgs.count_rows())
    tob_wgs_path = output_path('tob_wgs_filtered.mt')
    tob_wgs = tob_wgs.repartition(1000, shuffle=False)
    tob_wgs.write(tob_wgs_path)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
