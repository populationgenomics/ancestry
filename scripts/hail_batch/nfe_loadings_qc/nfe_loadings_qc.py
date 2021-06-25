"""Explore loadings dataset"""

import subprocess
from collections import Counter
import click
import hail as hl

LOADINGS = 'gs://cpg-tob-wgs-main/1kg_hgdp_densified_nfe/v0/loadings.ht/'
HGDP1KG_TOBWGS = (
    'gs://cpg-tob-wgs-main/1kg_hgdp_densified_pca/v2/'
    'hgdp1kg_tobwgs_joined_all_samples.mt/'
)


@click.command()
@click.option('--output', help='GCS output path', required=True)
def query(output):  # pylint: disable=too-many-locals
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    # get frequency of loadings values
    loadings = hl.read_table(LOADINGS)
    number_of_pcs = hl.len(loadings.loadings).take(1)[0]
    print(loadings.count())
    for i in range(0, (number_of_pcs)):
        pc = i + 1
        freq = Counter(hl.abs(loadings.loadings[i]).collect())
        filename = 'loadings_pc' + str(pc) + '.txt'
        with open(filename, 'w') as f:
            for key, value in freq.items():
                str_value = repr(key) + ' ' + repr(value)
                f.write(str_value + '\n')
        f.close()
        subprocess.run(['gsutil', 'cp', filename, output], check=False)

    # pull out variants that looked like they're capped in the loadings plot
    mt = hl.read_matrix_table(HGDP1KG_TOBWGS)
    # Get NFE samples only
    mt = mt.filter_cols(
        (mt.hgdp_1kg_metadata.population_inference.pop == 'nfe')
        | (mt.s.contains('TOB'))
    )
    intervals = [
        hl.parse_locus(x, reference_genome='GRCh38')
        for x in [
            'chr1:176163025',
            'chr5:272714',
            'chr5:36104012',
            'chr1:183565810',
            'chr3:58111799',
        ]
    ]
    mt_hits = mt.filter_rows(hl.literal(intervals).contains(mt.locus))
    mt_path = f'{output}/capped_loadings_intervals.mt'
    mt_hits.write(mt_path)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
