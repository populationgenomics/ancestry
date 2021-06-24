"""Explore loadings dataset"""

import subprocess
from collections import Counter
import click
import hail as hl

LOADINGS = 'gs://cpg-tob-wgs-test/1kg_hgdp_densify/v15/loadings.ht/'
# 'gs://cpg-tob-wgs-main/1kg_hgdp_densified_nfe/v0/loadings.ht/'


@click.command()
@click.option('--output', help='GCS output path', required=True)
def query(output):  # pylint: disable=too-many-locals
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    loadings = hl.read_table(LOADINGS)
    print(loadings.count())
    for i in range(0, 20):
        pc = i + 1
        freq = Counter(hl.abs(loadings.loadings[i]).collect())
        filename = 'loadings_pc' + str(pc) + '.txt'
        with open(filename, 'w') as f:
            for key, value in freq.items():
                str_value = repr(key) + ' ' + repr(value)
                f.write(str_value + '\n')
        f.close()
        subprocess.run(['gsutil', 'cp', filename, output], check=False)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
