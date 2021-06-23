"""Explore loadings dataset"""

from collections import Counter
import click
import hail as hl

LOADINGS = 'gs://cpg-tob-wgs-test/1kg_hgdp_densify/v15/loadings.ht/'
# 'gs://cpg-tob-wgs-main/1kg_hgdp_densified_nfe/v0/loadings.ht/'


@click.command()
# @click.option('--output', help='GCS output path', required=True)
def query():  # pylint: disable=too-many-locals
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    loadings = hl.read_table(LOADINGS)
    print(loadings.count())
    for i in range(0, 20):
        freq = Counter(hl.abs(loadings.loadings[i]).collect())
        print(freq)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
