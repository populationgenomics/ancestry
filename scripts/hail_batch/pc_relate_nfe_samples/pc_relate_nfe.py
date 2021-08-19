"""
Perform pc_relate on nfe samples from the HGDP/1KG dataset.
"""

import hail as hl


GNOMAD_HGDP_1KG_MT = (
    'gs://gcp-public-data--gnomad/release/3.1/mt/genomes/'
    'gnomad.genomes.v3.1.hgdp_1kg_subset_dense.mt'
)


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(GNOMAD_HGDP_1KG_MT)
    nrows_mt = mt.count_rows()
    mt = mt.sample_rows(100000 / nrows_mt, seed=12345)

    # Remove related samples (at the 2nd degree or closer)
    pc_rel = hl.pc_relate(mt.GT, 0.01, k=10, statistics='kin')
    known_trios = pc_rel.filter((pc_rel.i.s == 'HG01696') | (pc_rel.i.s == 'HG01629'))
    # Select parents from the child_pca_outlier matrix
    known_trios = known_trios.filter(
        (known_trios.j.s == 'HG01628')
        | (known_trios.j.s == 'HG01694')
        | (known_trios.j.s == 'HG01695')
        | (known_trios.j.s == 'HG01630')
    )
    print(known_trios.to_pandas())


if __name__ == '__main__':
    query()
