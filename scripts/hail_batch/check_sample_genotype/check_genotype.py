"""Get genotypes of samples"""

import click
import hail as hl

HGDP1KG_TOBWGS = (
    'gs://cpg-tob-wgs-analysis/1kg_hgdp_tobwgs_pca/'
    'v1/hgdp1kg_tobwgs_joined_all_samples.mt/'
)

GNOMAD_HGDP_1KG_MT = (
    'gs://gcp-public-data--gnomad/release/3.1/mt/genomes/'
    'gnomad.genomes.v3.1.hgdp_1kg_subset_dense.mt'
)


@click.command()
@click.option('--output', help='GCS output path', required=True)
def query(output):  # pylint: disable=too-many-locals
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    # compute allele frequency count on amr population
    mt = hl.read_matrix_table(HGDP1KG_TOBWGS)
    amr = mt.filter_cols(mt.hgdp_1kg_metadata.population_inference.pop == 'amr')
    amr = amr.annotate_rows(gt_stats=hl.agg.call_stats(amr.GT, amr.alleles))
    amr_af_path = f'{output}/amr.tsv'
    amr.gt_stats.AF.export(amr_af_path)

    # get TOB-WGS GT data for downstream manipulation
    mt_gt_path = f'{output}/tob_wgs_GT.tsv'
    mt.GT.export(mt_gt_path)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
