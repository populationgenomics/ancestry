"""
Test densify function on TOB-WGS data.
"""

import click
import hail as hl
from hail.experimental import lgt_to_gt


GNOMAD_LIFTOVER_LOADINGS = 'gs://cpg-reference/gnomad/gnomad_loadings_90k_liftover.ht'

GNOMAD_HGDP_1KG_MT = (
    'gs://gcp-public-data--gnomad/release/3.1/mt/genomes/'
    'gnomad.genomes.v3.1.hgdp_1kg_subset_dense.mt'
)

TOB_WGS = 'gs://cpg-tob-wgs-main/joint_vcf/v1/raw/genomes.mt'


@click.command()
@click.option('--output', help='GCS output path', required=True)
def query(output):  # pylint: disable=too-many-locals
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    # get TOB-WGS allele frequencies
    tob_wgs = hl.read_matrix_table(TOB_WGS).key_rows_by('locus', 'alleles')
    tob_wgs = hl.experimental.densify(tob_wgs)
    tob_wgs = tob_wgs.annotate_entries(GT=lgt_to_gt(tob_wgs.LGT, tob_wgs.LA))
    tob_wgs = tob_wgs.annotate_rows(
        gt_stats=hl.agg.call_stats(tob_wgs.GT, tob_wgs.alleles)
    )

    # Get gnomAD allele frequency of variants that aren't in TOB-WGS
    loadings_gnomad = hl.read_table(GNOMAD_LIFTOVER_LOADINGS).key_by('locus', 'alleles')
    hgdp_1kg = hl.read_matrix_table(GNOMAD_HGDP_1KG_MT)
    hgdp_1kg_row = hgdp_1kg.rows()[loadings_gnomad.locus, loadings_gnomad.alleles]
    tob_wgs_row = tob_wgs.rows()[loadings_gnomad.locus, loadings_gnomad.alleles]
    loadings_gnomad = loadings_gnomad.annotate(
        gnomad_AF=hgdp_1kg_row.gnomad_freq.AF,
        gnomad_popmax_AF=hgdp_1kg_row.gnomad_popmax.AF,
        TOB_WGS_AF=tob_wgs_row.gt_stats.AF,
    )
    population_af_metadata = hgdp_1kg.gnomad_freq_meta.collect()
    loadings_gnomad = loadings_gnomad.annotate_globals(
        gnomad_freq_meta=population_af_metadata
    )
    gnomad_variants = loadings_gnomad.drop('loadings')
    gnomad_variants_path = f'{output}/gnomad_annotated_variants_densified.ht'
    gnomad_variants.write(gnomad_variants_path)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
