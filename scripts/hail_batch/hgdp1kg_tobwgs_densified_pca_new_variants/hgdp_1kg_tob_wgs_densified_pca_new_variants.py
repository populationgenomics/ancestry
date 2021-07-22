"""
Perform PCA on densified TOB-WGS data. Reliant on output from
```
hgdp1kg_tobwgs_densified_pca_new_variants/
hgdp_1kg_tob_wgs_densified_pca_new_variants.py
````
"""

import hail as hl
import pandas as pd
from hail.experimental import lgt_to_gt
from analysis_runner import bucket_path, output_path


# TOB_WGS = bucket_path('1kg_hgdp_densify_new_variants/v0/tob_wgs_filtered.mt/')
TOB_WGS = bucket_path('mt/v4.mt/')
# TOB_WGS = 'gs://cpg-tob-wgs-test/mt/v4.mt/'
GNOMAD_HGDP_1KG_MT = (
    'gs://gcp-public-data--gnomad/release/3.1/mt/genomes/'
    'gnomad.genomes.v3.1.hgdp_1kg_subset_dense.mt'
)


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    tob_wgs = hl.read_matrix_table(TOB_WGS).head(10000)
    hgdp_1kg = hl.read_matrix_table(GNOMAD_HGDP_1KG_MT).head(10000)

    # keep loci that are contained in the densified, filtered tob-wgs mt
    hgdp_1kg = hgdp_1kg.semi_join_rows(tob_wgs.rows())

    # Entries and columns must be identical
    tob_wgs_select = tob_wgs.select_entries(
        GT=lgt_to_gt(tob_wgs.LGT, tob_wgs.LA)
    ).select_cols()
    hgdp_1kg_select = hgdp_1kg.select_entries(hgdp_1kg.GT).select_cols()
    # Join datasets
    hgdp1kg_tobwgs_joined = hgdp_1kg_select.union_cols(tob_wgs_select)
    # Add in metadata information
    hgdp_1kg_metadata = hgdp_1kg.cols()
    hgdp1kg_tobwgs_joined = hgdp1kg_tobwgs_joined.annotate_cols(
        hgdp_1kg_metadata=hgdp_1kg_metadata[hgdp1kg_tobwgs_joined.s]
    )
    # save this for population-level PCAs
    mt_path = output_path('hgdp1kg_tobwgs_joined_all_samples.mt')
    if not hl.hadoop_exists(mt_path):
        hgdp1kg_tobwgs_joined.write(mt_path)

    # Perform PCA
    eigenvalues_path = output_path('eigenvalues.ht')
    scores_path = output_path('scores.ht')
    loadings_path = output_path('loadings.ht')
    eigenvalues, scores, loadings = hl.hwe_normalized_pca(
        hgdp1kg_tobwgs_joined.GT, compute_loadings=True, k=20
    )
    hl.Table.from_pandas(pd.DataFrame(eigenvalues)).export(eigenvalues_path)
    scores.write(scores_path, overwrite=True)
    loadings.write(loadings_path, overwrite=True)


if __name__ == '__main__':
    query()
