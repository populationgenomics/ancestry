"""Get genotypes of samples"""

import click
import pandas as pd
from pyarrow import feather
import hail as hl

HGDP1KG_TOBWGS = (
    'gs://cpg-tob-wgs-analysis/1kg_hgdp_tobwgs_pca/v0/hgdp1kg_tobwgs_joined.mt/'
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
    hgdp_1kg = hl.read_matrix_table(GNOMAD_HGDP_1KG_MT)
    amr = hgdp_1kg.filter_cols(hgdp_1kg.population_inference.pop == 'amr')
    amr = amr.annotate_rows(gt_stats=hl.agg.call_stats(amr.GT, amr.alleles))
    amr_to_pandas = amr.rows().to_pandas()
    amr_af = pd.DataFrame(amr_to_pandas['gt_stats.AF'])
    amr_af_path = f'{output}/amr_AF.feather'
    feather.write_feather(
        amr_af,
        amr_af_path,
    )
    # get three TOB-WGS samples
    mt = hl.read_matrix_table(HGDP1KG_TOBWGS)
    mt_filt = mt.filter_cols(
        (mt.s == 'TOB1524') | (mt.s == 'TOB1532') | (mt.s == 'TOB1533')
    )
    mt_filt = hl.sample_qc(mt_filt)
    mt_filt_to_pandas = mt_filt.cols().to_pandas()
    mt_filt_path = f'{output}/tob_wgs_selected_sample_qc.csv'
    mt_filt_to_pandas.to_csv(mt_filt_path, index=False)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
