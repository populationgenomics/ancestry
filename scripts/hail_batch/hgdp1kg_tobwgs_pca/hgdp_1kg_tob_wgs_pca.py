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
        hl.is_defined(loadings.index(hgdp_1kg['locus'], hgdp_1kg['alleles'])) &
        hl.is_defined(tob_wgs.index(hgdp_1kg['locus'], hgdp_1kg['alleles']))
    )
    tob_wgs = tob_wgs.semi_join_rows(hgdp_1kg)

    # Join datasets by merging columns
    # Entries and columns must be identical
    change_lgt_to_gt = join_rows_tobwgs_liftover.annotate_entries(
        GT=lgt_to_gt(join_rows_tobwgs_liftover.LGT, join_rows_tobwgs_liftover.LA)
    )
    select_entries_hgdp1kg = join_rows_hgdp1kg_liftover.select_entries(
        join_rows_hgdp1kg_liftover.GT
    )
    select_entries_tobwgs = change_lgt_to_gt.select_entries(change_lgt_to_gt.GT)

    # Only keep columns we need
    select_columns_hgdp1kg = select_entries_hgdp1kg.select_cols(
        select_entries_hgdp1kg.sample_qc
    )
    hgdp1kg_filtered = select_columns_hgdp1kg.drop('sample_qc')
    tobwgs_filtered = select_entries_tobwgs
    # Join datasets
    hgdp1kg_tobwgs_joined = hgdp1kg_filtered.union_cols(
        tobwgs_filtered, row_join_type='inner'
    )

    # Perform PCA
    eigenvalues_path = f'{output}/eigenvalues.csv'
    scores_path = f'{output}/scores.ht'
    loadings_path = f'{output}/loadings.ht'

    # test on 100 samples
    mt_head = hgdp1kg_tobwgs_joined.head(
        n=hgdp1kg_tobwgs_joined.count_rows(), n_cols=100
    )
    eigenvalues, scores, loadings = hl.hwe_normalized_pca(
        mt_head.GT, compute_loadings=True, k=20
    )
    # save the list of eigenvalues
    eigenvalues_df = pd.DataFrame(eigenvalues)
    eigenvalues_df.to_csv(eigenvalues_path, index=False)
    # save the scores and loadings as a hail table
    scores.write(scores_path, overwrite=True)
    loadings.write(loadings_path, overwrite=True)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
