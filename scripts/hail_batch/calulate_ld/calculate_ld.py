"""Calculate ld using the ld_matrix function"""

import hail as hl

TOB_WGS = "gs://cpg-tob-wgs-test/mt/v7.mt/"


def query():
    """Query script entry point."""

    hl.init(default_reference="GRCh38")

    tob_wgs = hl.read_matrix_table(TOB_WGS)
    tob_wgs = hl.experimental.densify(tob_wgs)
    print("Starting ld calculation")
    ld = hl.ld_matrix(tob_wgs.GT.n_alt_alleles(), mt.locus, radius=1e6)
    ld.to_numpy()


if __name__ == "__main__":
    query()
