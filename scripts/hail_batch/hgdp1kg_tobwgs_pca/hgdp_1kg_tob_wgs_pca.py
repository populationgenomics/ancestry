"""Generates loadings, scores, and eigenvalues for the HGDP,1KG, and tob-wgs dataset"""

import click
import pandas as pd
import hail as hl
from hail.experimental import lgt_to_gt

GNOMAD_HGDP_1KG_MT = (
    'gs://gcp-public-data--gnomad/release/3.1/mt/genomes/'
    'gnomad.genomes.v3.1.hgdp_1kg_subset_dense.mt'
)

TOB_WGS = 'gs://cpg-tob-wgs-temporary/joint_vcf/v1/work/combiner/genomes.mt'

GNOMAD_V2_LOADINGS = (
    'gs://gcp-public-data--gnomad/release/2.1/pca/gnomad.r2.1.pca_loadings.ht'
)


@click.command()
@click.option('--output', help='GCS output path', required=True)
def query(output):  # pylint: disable=too-many-locals
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    # Filter the matrix tables to rows whose keys appear in both datasets
    mt_hgdp_1kg = hl.read_matrix_table(GNOMAD_HGDP_1KG_MT)
    mt_tob_wgs = hl.read_matrix_table(TOB_WGS)
    # Remove 'allele' row key so that rows are the same (necessary step to join)
    mt_hgdp_1kg = mt_hgdp_1kg.key_rows_by('locus')
    mt_hgdp_1kg = mt_hgdp_1kg.drop('alleles')
    join_rows_hgdp1kg_tobwgs = mt_hgdp_1kg.semi_join_rows(mt_tob_wgs.rows())
    join_rows_tobwgs_hgdp1kg = mt_tob_wgs.semi_join_rows(mt_hgdp_1kg.rows())

    # Join dataset by merging columns
    # Entries and columns must be identical
    change_lgt_to_gt = join_rows_tobwgs_hgdp1kg.annotate_entries(
        GT=lgt_to_gt(join_rows_tobwgs_hgdp1kg.LGT, join_rows_tobwgs_hgdp1kg.LA)
    )
    select_entries_hgdp1kg = join_rows_hgdp1kg_tobwgs.select_entries(
        join_rows_hgdp1kg_tobwgs.GT
    )
    select_entries_tobwgs = change_lgt_to_gt.select_entries(change_lgt_to_gt.GT)
    # Only keep columns we need
    select_columns_hgdp1kg = select_entries_hgdp1kg.select_cols(
        select_entries_hgdp1kg.sample_qc
    )
    hgdp1kg_filtered = select_columns_hgdp1kg.drop('sample_qc')
    tobwgs_filtered = select_entries_tobwgs
    # Join datasets
    hgdp1kg_tobwgs_joined = hgdp1kg_filtered.union_cols(tobwgs_filtered)

    # liftover and get variants
    ht_gnomad_loadings = hl.read_table(GNOMAD_V2_LOADINGS)
    ht_gnomad_loadings.count()
    # 94176
    rg37 = hl.get_reference('GRCh37')
    rg38 = hl.get_reference('GRCh38')
    rg37.add_liftover(
        'gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38
    )
    ht_gnomad_loadings_liftover = ht_gnomad_loadings.annotate(
        liftover=hl.liftover(ht_gnomad_loadings.locus, 'GRCh38', include_strand=False),
        old_locus=ht_gnomad_loadings.locus,
    )
    ht_gnomad_loadings_liftover = ht_gnomad_loadings_liftover.key_by(
        locus=ht_gnomad_loadings_liftover.liftover
    )

    # add liftover into the filtered, combined matrix table
    hgdp1kg_tobwgs_liftover = hgdp1kg_tobwgs_joined.annotate_rows(
        rows_to_keep=ht_gnomad_loadings_liftover[hgdp1kg_tobwgs_joined.locus]
    )
    hgdp1kg_tobwgs_filt = hgdp1kg_tobwgs_liftover.filter_rows(
        hgdp1kg_tobwgs_liftover.locus == hgdp1kg_tobwgs_liftover.rows_to_keep.liftover
    )

    # Perform PCA
    eigenvalues_path = f'{output}/eigenvalues.csv'
    scores_path = f'{output}/scores.ht'
    loadings_path = f'{output}/loadings.ht'

    # test on 100 samples
    mt_head = hgdp1kg_tobwgs_filt.head(n=hgdp1kg_tobwgs_filt.count_rows(), n_cols=100)
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
