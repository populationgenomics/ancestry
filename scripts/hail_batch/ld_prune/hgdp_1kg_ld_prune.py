"""Generates ld pruning for the HGDP + 1KG dataset"""

import click
from gnomad.utils.annotations import annotate_adj
import hail as hl

GNOMAD_HGDP_1KG_MT = (
    'gs://gcp-public-data--gnomad/release/3.1/mt/genomes/'
    'gnomad.genomes.v3.1.hgdp_1kg_subset_dense.mt'
)


@click.command()
@click.option('--output', help='GCS output path', required=True)
def query(output):
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    mt_path = f'{output}/filtered_mt.mt'
    mt = hl.read_matrix_table(GNOMAD_HGDP_1KG_MT)
    # reproduce gnomAD genotype filtering
    mt = annotate_adj(mt)
    mt = mt.filter_entries(mt.adj)
    # perform ld pruning
    biallelic_mt = mt.filter_rows(hl.len(mt.alleles) == 2)
    pruned_variant_table = hl.ld_prune(biallelic_mt.GT, r2=0.2, bp_window_size=500000)
    filtered_mt = mt.filter_rows(hl.is_defined(pruned_variant_table[mt.row_key]))
    # save filtered mt table
    filtered_mt.write(mt_path, overwrite=True)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
