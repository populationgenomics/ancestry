"""Tests the densify function on the tob-wgs dataset"""
# As a first step, could we just do the densify before filter (that seems safer),

import click
import hail as hl

GNOMAD_HGDP_1KG_MT = (
    'gs://gcp-public-data--gnomad/release/3.1/mt/genomes/'
    'gnomad.genomes.v3.1.hgdp_1kg_subset_dense.mt'
)

TOB_WGS = 'gs://cpg-tob-wgs-test/mt/test-v1-raw.mt'

GNOMAD_LIFTOVER_LOADINGS = 'gs://cpg-reference/gnomad/gnomad_loadings_90k_liftover.ht'


@click.command()
@click.option('--output', help='GCS output path', required=True)
def query(output):  # pylint: disable=too-many-locals
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    hgdp_1kg = hl.read_matrix_table(GNOMAD_HGDP_1KG_MT)
    tob_wgs = hl.read_matrix_table(TOB_WGS).key_rows_by('locus', 'alleles')
    loadings = hl.read_table(GNOMAD_LIFTOVER_LOADINGS).key_by('locus', 'alleles')

    # filter to loci that are contained in both tables and the loadings after densifying
    tob_wgs = hl.experimental.densify(tob_wgs)
    hgdp_1kg = hgdp_1kg.filter_rows(
        hl.is_defined(loadings.index(hgdp_1kg['locus'], hgdp_1kg['alleles']))
        & hl.is_defined(tob_wgs.index_rows(hgdp_1kg['locus'], hgdp_1kg['alleles']))
    )
    tob_wgs = tob_wgs.semi_join_rows(hgdp_1kg.rows())
    tob_wgs_path = f'{output}/tob_wgs_filtered_densified.mt'
    tob_wgs = tob_wgs.checkpoint(tob_wgs_path)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
