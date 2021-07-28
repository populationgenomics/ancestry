"""
Perform pca on samples specific to a population
from the HGDP,1KG, and tob-wgs dataset after densifying.
Reliant on output from

```hgdp1kg_tobwgs_densified_pca_new_variants/
hgdp_1kg_tob_wgs_densified_pca_new_variants.py
```
"""

import hail as hl
import click


HGDP1KG_TOBWGS = (
    'gs://cpg-tob-wgs-test/1kg_hgdp_densified_pca_new_variants/'
    'v0/hgdp1kg_tobwgs_joined_all_samples.mt'
)


@click.command()
@click.option('--pop', help='Population to subset from the 1KG (e.g. afr, nfe)')
def query(pop):
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(HGDP1KG_TOBWGS)
    if pop:
        # Get samples from the specified population only
        mt = mt.filter_cols(
            (mt.hgdp_1kg_metadata.population_inference.pop == pop.lower())
            | (mt.s.contains('TOB'))
        )
    else:
        mt = mt.filter_cols(mt.s.contains('TOB'))

    # Perform PCA
    eigenvalues, scores, loadings = hl.hwe_normalized_pca(
        mt.GT, compute_loadings=True, k=20
    )
    print(eigenvalues)
    print(scores)
    print(loadings)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
