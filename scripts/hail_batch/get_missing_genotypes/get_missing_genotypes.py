"""Get number of variants with a MAF > 0.05"""

import hail as hl
from analysis_runner import bucket_path

TOB_WGS = bucket_path('mt/v5.1.mt/')


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    tob_wgs = hl.read_matrix_table(TOB_WGS)
    # densify the table and split multiallelics
    tob_wgs = hl.experimental.densify(tob_wgs)
    tob_wgs = hl.split_multi_hts(tob_wgs)
    # only keep biallelic and multiallelic loci
    tob_wgs = tob_wgs.filter_rows(hl.len(tob_wgs.alleles) >= 2)
    # count number of rows before filtering and print
    nrows_pre = tob_wgs.count_rows()
    print(f'There are {nrows_pre} rows before filtering')
    # filter out rows with NA genotypes
    tob_wgs = tob_wgs.filter_rows(hl.agg.any(hl.is_missing(tob_wgs.GT)))
    nrows_post = tob_wgs.count_rows()
    print(f'There are {nrows_post} rows after filtering')


if __name__ == '__main__':
    query()
