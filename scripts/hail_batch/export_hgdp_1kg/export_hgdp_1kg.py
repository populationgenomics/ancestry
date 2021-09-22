"""Export HGDP/1kG subset as PLINK format"""

import hail as hl
from analysis_runner import output_path

GNOMAD_HGDP_1KG_MT = (
    'gs://gcp-public-data--gnomad/release/3.1/mt/genomes/'
    'gnomad.genomes.v3.1.hgdp_1kg_subset_dense.mt'
)


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(GNOMAD_HGDP_1KG_MT)
    samples = [
        'HG00100',
        'NA20773',
        'HG00102',
        'HG01513',
        'v3.1::NA12878',
        'v3.1::NA12891',
        'NA12892',
    ]
    mt = mt.filter_cols(hl.literal(samples).contains(mt.s), keep=True)
    nrows = mt.count_rows()
    mt = mt.sample_rows(10000 / nrows, seed=12345)
    mt = hl.split_multi_hts(mt)
    mt_path = output_path('hgdp_1kg_plink_subset')
    hl.export_plink(mt, mt_path, ind_id=mt.s)


if __name__ == '__main__':
    query()
