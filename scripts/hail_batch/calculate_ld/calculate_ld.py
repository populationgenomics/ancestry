"""Calculate ld using the ld_matrix function"""

import hail as hl
import pandas as pd
from analysis_runner import bucket_path, output_path

TOB_WGS = bucket_path('mt/v7.mt/')


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    tob_wgs = hl.read_matrix_table(TOB_WGS)
    tob_wgs = hl.experimental.densify(tob_wgs)
    # filter out constant variants
    tob_wgs = tob_wgs.filter_rows(hl.len(tob_wgs.alleles) == 2)
    tob_wgs = tob_wgs.head(30000)
    ld = hl.ld_matrix(tob_wgs.GT.n_alt_alleles(), tob_wgs.locus, radius=1e6)
    ld = pd.DataFrame(ld.to_numpy())
    # save pandas df
    ld_filename = output_path(f'ld_matrix.csv', 'analysis')
    ld.to_csv(ld_filename, index=False)


if __name__ == '__main__':
    query()
