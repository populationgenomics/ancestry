"""
Perform pc_relate on nfe samples from the HGDP/1KG dataset.
"""

import hail as hl


GNOMAD_HGDP_1KG_MT = (
    'gs://gcp-public-data--gnomad/release/3.1/mt/genomes/'
    'gnomad.genomes.v3.1.hgdp_1kg_subset_dense.mt'
)

GNOMAD_V2_LOADINGS = (
    'gs://gcp-public-data--gnomad/release/2.1/pca/gnomad.r2.1.pca_loadings.ht'
)


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(GNOMAD_HGDP_1KG_MT)
    # Get samples from the specified population only
    mt = mt.filter_cols(mt.population_inference.pop == 'nfe')
    mt = mt.head(1000)

    # liftover and get variants
    ht_gnomad_loadings = hl.read_table(GNOMAD_V2_LOADINGS)
    rg37 = hl.get_reference('GRCh37')
    rg38 = hl.get_reference('GRCh38')
    rg37.add_liftover(
        'gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38
    )
    ht_gnomad_loadings_liftover = ht_gnomad_loadings.annotate(
        liftover=hl.liftover(ht_gnomad_loadings.locus, 'GRCh38', include_strand=False),
        old_locus=ht_gnomad_loadings.locus,
    )
    ht_gnomad_loadings_liftover = ht_gnomad_loadings_liftover.key_by(
        locus=ht_gnomad_loadings_liftover.liftover,
        alleles=ht_gnomad_loadings_liftover.alleles,
    )
    mt = mt.semi_join_rows(ht_gnomad_loadings_liftover)

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
