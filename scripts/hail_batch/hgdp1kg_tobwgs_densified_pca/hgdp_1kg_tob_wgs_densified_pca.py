"""
Perform PCA on densified TOB-WGS data.
Reliant on output from hgdp_1kg_tob_wgs_densify.py
"""

import click
import pandas as pd
import hail as hl
from hail.experimental import lgt_to_gt

# contains batches 1-4
TOB_WGS = 'gs://cpg-tob-wgs-main/1kg_hgdp_densify/v5/tob_wgs_filtered.mt/'

GNOMAD_HGDP_1KG_MT = (
    'gs://gcp-public-data--gnomad/release/3.1/mt/genomes/'
    'gnomad.genomes.v3.1.hgdp_1kg_subset_dense.mt'
)


@click.command()
@click.option('--output', help='GCS output path', required=True)
def query(output):  # pylint: disable=too-many-locals
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    tob_wgs = hl.read_matrix_table(TOB_WGS)
    hgdp_1kg = hl.read_matrix_table(GNOMAD_HGDP_1KG_MT)

    # keep loci that are contained in the densified, filtered tob-wgs mt
    hgdp_1kg = hgdp_1kg.semi_join_rows(tob_wgs.rows())

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
    eigenvalues_path = f'{output}/eigenvalues.ht'
    scores_path = f'{output}/scores.ht'
    loadings_path = f'{output}/loadings.ht'
    eigenvalues, scores, loadings = hl.hwe_normalized_pca(
        hgdp1kg_tobwgs_joined.GT, compute_loadings=True, k=20
    )
    # save the list of eigenvalues
    hl.Table.from_pandas(pd.DataFrame(eigenvalues)).export(eigenvalues_path)
    # save the scores and loadings as a hail table
    scores.write(scores_path, overwrite=True)
    loadings.write(loadings_path, overwrite=True)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
