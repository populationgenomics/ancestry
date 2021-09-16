"""Repartition lifted over gnomAD v2 loadings"""

import hail as hl
from analysis_runner import output_path

LOADINGS = 'gs://cpg-reference/gnomad/gnomad_loadings_90k_liftover.ht'


def query():
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    loadings = hl.read_table(LOADINGS)
    loadings = loadings.repartition(100, shuffle=False)
    loadings_path = output_path(f'gnomad_loadings_90k_liftover_repartitioned.ht')
    loadings.write(loadings_path)


if __name__ == '__main__':
    query()
