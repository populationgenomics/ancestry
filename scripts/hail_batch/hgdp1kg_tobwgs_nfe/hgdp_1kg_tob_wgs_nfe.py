"""Perform pca on nfe samples from the HGDP,1KG, and tob-wgs dataset"""

import click
import pandas as pd
import hail as hl
from hail.experimental import lgt_to_gt

HGDP1KG_TOBWGS = (
    'gs://cpg-tob-wgs-analysis/1kg_hgdp_tobwgs_pca/'
    'v1/hgdp1kg_tobwgs_joined_all_samples.mt'
)

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

    mt = hl.read_matrix_table(HGDP1KG_TOBWGS)
    # Get NFE samples only
    mt = mt.filter_cols(
        (mt.hgdp_1kg_metadata.population_inference.pop == 'nfe')
        | (mt.s.contains('TOB'))
    )

    # Perform PCA
    eigenvalues_path = f'{output}/eigenvalues.csv'
    scores_path = f'{output}/scores.ht'
    loadings_path = f'{output}/loadings.ht'
    eigenvalues, scores, loadings = hl.hwe_normalized_pca(
        mt.GT, compute_loadings=True, k=20
    )
    eigenvalues_df = pd.DataFrame(eigenvalues)
    eigenvalues_df.to_csv(eigenvalues_path, index=False)
    scores.write(scores_path, overwrite=True)
    loadings.write(loadings_path, overwrite=True)

    # get TOB-WGS allele frequencies
    tob_wgs = hl.read_matrix_table(TOB_WGS).key_rows_by('locus', 'alleles')
    tob_wgs = tob_wgs.annotate_entries(GT=lgt_to_gt(tob_wgs.LGT, tob_wgs.LA))
    tob_wgs = tob_wgs.annotate_rows(
        gt_stats=hl.agg.call_stats(tob_wgs.GT, tob_wgs.alleles)
    )

    # Get gnomAD allele frequency of variants that aren't in TOB-WGS
    loadings_gnomad = hl.read_table(GNOMAD_LIFTOVER_LOADINGS).key_by('locus', 'alleles')
    hgdp_1kg = hl.read_matrix_table(GNOMAD_HGDP_1KG_MT)
    loadings_gnomad = loadings_gnomad.annotate(
        gnomad_AF=hgdp_1kg.rows()[
            loadings_gnomad.locus, loadings_gnomad.alleles
        ].gnomad_freq.AF,
        gnomad_popmax_AF=hgdp_1kg.rows()[
            loadings_gnomad.locus, loadings_gnomad.alleles
        ].gnomad_popmax.AF,
        TOB_variant=hl.is_defined(
            mt.rows()[loadings_gnomad.locus, loadings_gnomad.alleles].gnomad_popmax.AF
        ),
        TOB_WGS_AF=tob_wgs.rows()[
            loadings_gnomad.locus, loadings_gnomad.alleles
        ].gt_stats.AF,
    )
    loadings_gnomad_path = f'{output}/gnomad_loadings_annotated_variants.mt'
    mt.write(loadings_gnomad_path)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
