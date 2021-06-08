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

# contains batches 1-4
TOB_WGS = 'gs://cpg-tob-wgs-test/mt/test-v1-raw.mt'


@click.command()
@click.option('--output', help='GCS output path', required=True)
def query(output):  # pylint: disable=too-many-locals
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    hgdp_1kg = hl.read_matrix_table(GNOMAD_HGDP_1KG_MT)
    tob_wgs = hl.read_matrix_table(TOB_WGS).key_rows_by('locus', 'alleles')
    loadings = hl.read_table(GNOMAD_LIFTOVER_LOADINGS).key_by('locus', 'alleles')

    # filter to loci that are contained in both tables and the loadings after densifying
    hgdp_1kg = hgdp_1kg.filter_rows(hgdp_1kg.locus.contig == 'chr22')
    hgdp_1kg = hgdp_1kg.head(1000000)
    tob_wgs = tob_wgs.filter_rows(tob_wgs.locus.contig == 'chr22')
    tob_wgs = tob_wgs.head(1000000)
    tob_wgs = hl.experimental.densify(tob_wgs)
    hgdp_1kg = hgdp_1kg.filter_rows(
        hl.is_defined(loadings.index(hgdp_1kg['locus'], hgdp_1kg['alleles']))
        & hl.is_defined(tob_wgs.index_rows(hgdp_1kg['locus'], hgdp_1kg['alleles']))
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

    # Perform PCA
    scores_path = f'{output}/scores.ht'
    loadings_path = f'{output}/loadings.ht'
    _, scores, loadings = hl.hwe_normalized_pca(
        hgdp1kg_tobwgs_joined.GT, compute_loadings=True, k=20
    )
    # save the scores and loadings as a hail table
    scores.write(scores_path, overwrite=True)
    loadings.write(loadings_path, overwrite=True)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
