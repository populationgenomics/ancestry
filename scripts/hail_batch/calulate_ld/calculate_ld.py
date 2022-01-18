"""Calculate ld using the ld_matrix function"""

import hail as hl
import pandas as pd

TOB_WGS = 'gs://cpg-tob-wgs-test/mt/v7.mt/'


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    tob_wgs = hl.read_matrix_table(TOB_WGS)
    tob_wgs = hl.experimental.densify(tob_wgs)
    # filter out constant variants
    tob_wgs = tob_wgs.filter_rows(hl.len(tob_wgs.alleles) == 2)
    tob_wgs = tob_wgs.head(40000)
    ld = hl.ld_matrix(tob_wgs.GT.n_alt_alleles(), tob_wgs.locus, radius=2e6)
    ld = ld.to_numpy()
    ld = pd.DataFrame(ld)
    # save pandas df
    ld.to_csv('gs://cpg-tob-wgs-test/kat/v1/ld_matrix_40k.csv')
    # get number of non-zero, non na values across rows
    nonzero_values = ld.fillna(0).astype(bool).sum(axis=1)
    nonzero_values.to_csv(
        'gs://cpg-tob-wgs-test/kat/v1/nonzero_nona_values_2M_radius_40k.csv'
    )
    # get number of positive values (including 1's) across rows
    positive_values = ld.fillna(0).gt(0).sum(axis=1)
    positive_values.to_csv(
        'gs://cpg-tob-wgs-test/kat/v1/positive_values_2M_radius_40k.csv'
    )
    # get number of negative values across rows
    negative_values = ld.fillna(0).lt(0).sum(axis=1)
    negative_values.to_csv(
        'gs://cpg-tob-wgs-test/kat/v1/negative_values_2M_radius_40k.csv'
    )


if __name__ == '__main__':
    query()
