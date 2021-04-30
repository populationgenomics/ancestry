"""Generates loadings, scores, and eigenvalues for the HGDP,1KG, and tob-wgs dataset"""

import click
import pandas as pd
import hail as hl
from hail.experimental import lgt_to_gt

GNOMAD_HGDP_1KG_MT = (
    'gs://gcp-public-data--gnomad/release/3.1/mt/genomes/'
    'gnomad.genomes.v3.1.hgdp_1kg_subset_dense.mt'
)

TOB_WGS = 'gs://cpg-tob-wgs-main/joint_vcf/v1/raw/genomes.mt'

GNOMAD_LIFTOVER_LOADINGS = 'gs://cpg-reference/gnomad/gnomad_loadings_90k_liftover.ht'


@click.command()
@click.option('--output', help='GCS output path', required=True)
def query(output):  # pylint: disable=too-many-locals
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    hgdp_1kg = hl.read_matrix_table(GNOMAD_HGDP_1KG_MT)
    tob_wgs = hl.read_matrix_table(TOB_WGS).key_rows_by('locus', 'alleles')
    loadings = hl.read_table(GNOMAD_LIFTOVER_LOADINGS).key_by('locus', 'alleles')

    # filter to loci that are contained in both tables and the loadings
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
    mt_path = f'{output}/hgdp1kg_tobwgs_joined_all_samples.mt'
    if not hl.hadoop_exists(mt_path):
        hgdp1kg_tobwgs_joined.write(mt_path)
    hgdp1kg_tobwgs_joined = hl.read_matrix_table(mt_path)

    # Perform PCA
    eigenvalues_path = f'{output}/eigenvalues.csv'
    scores_path = f'{output}/scores.ht'
    loadings_path = f'{output}/loadings.ht'
    eigenvalues, scores, loadings = hl.hwe_normalized_pca(
        hgdp1kg_tobwgs_joined.GT, compute_loadings=True, k=20
    )
    # save the list of eigenvalues
    eigenvalues_df = pd.DataFrame(eigenvalues)
    eigenvalues_df.to_csv(eigenvalues_path, index=False)
    # save the scores and loadings as a hail table
    scores.write(scores_path, overwrite=True)
    loadings.write(loadings_path, overwrite=True)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
