"""
Test densify function on TOB-WGS data.
"""

import click
import hail as hl


GNOMAD_LIFTOVER_LOADINGS = 'gs://cpg-reference/gnomad/gnomad_loadings_90k_liftover.ht'

GNOMAD_HGDP_1KG_MT = (
    'gs://gcp-public-data--gnomad/release/3.1/mt/genomes/'
    'gnomad.genomes.v3.1.hgdp_1kg_subset_dense.mt'
)

# contains batches 1-4
TOB_WGS = 'gs://cpg-tob-wgs-main/mt/v2-raw.mt/'


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
    tob_wgs_path = f'{output}/tob_wgs_filtered.mt'
    tob_wgs = tob_wgs.repartition(10, shuffle=False)
    tob_wgs.write(tob_wgs_path)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
