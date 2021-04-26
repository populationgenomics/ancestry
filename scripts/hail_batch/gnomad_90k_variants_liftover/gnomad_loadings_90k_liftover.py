"""Perform liftover for the gnomAD 90k variants"""

import click
import hail as hl

GNOMAD_V2_LOADINGS = (
    'gs://gcp-public-data--gnomad/release/2.1/pca/gnomad.r2.1.pca_loadings.ht'
)


@click.command()
@click.option('--output', help='GCS output path', required=True)
def query(output):  # pylint: disable=too-many-locals
    """Query script entry point."""

    hl.init(default_reference='GRCh38')

    gnomad_loadings_path = f'{output}/gnomad_loadings_90k_liftover.ht'

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
        locus=ht_gnomad_loadings_liftover.liftover
    )

    # save gnomad loadings
    ht_gnomad_loadings_liftover.write(gnomad_loadings_path, overwrite=True)


if __name__ == '__main__':
    query()  # pylint: disable=no-value-for-parameter
