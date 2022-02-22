#!/usr/bin/env python3

"""Entry point for the analysis runner."""

import os
import click
import hailtop.batch as hb
from analysis_runner import dataproc


@click.command()
@click.option('--script', 'script', help='path to main script')
@click.option(
    '--expression', required=True, help='A sample x gene TSV of expression values'
)
@click.option('--genotype', required=True, help='A TSV of genotypes for each sample')
@click.option(
    '--geneloc', required=True, help='A TSV of start and end positions for each gene'
)
@click.option(
    '--snploc',
    required=True,
    help='A TSV of snp IDs with chromsome and position values for each',
)
@click.option(
    '--covariates', required=True, help='A TSV of covariates to calculate residuals'
)
@click.option(
    '--keys',
    required=True,
    help='A TSV of sample ids to convert external to internal IDs',
)  # pylint: disable=too-many-locals
@click.option(
    '--output-prefix',
    required=True,
    help='A path prefix of where to output files, eg: gs://MyBucket/output-folder/',
)
def main(
    script: str,
    expression: str,
    genotype,
    geneloc,
    snploc,
    covariates,
    keys,
    output_prefix: str,
):
    """
    runs a script inside dataproc to execute generate_eqtl_spearman.py
    """

    service_backend = hb.ServiceBackend(
        billing_project=os.getenv('HAIL_BILLING_PROJECT'),
        bucket=os.getenv('HAIL_BUCKET'),
    )

    batch = hb.Batch(name='run-eqtl-association', backend=service_backend)

    job = dataproc.hail_dataproc_job(
        batch=batch,
        max_age='12h',
        num_secondary_workers=20,
        init=['gs://cpg-reference/hail_dataproc/install_common.sh'],
        job_name=f'run-eqtl-association',
        script=f'{script} --expression {expression} --genotype {genotype} --geneloc {geneloc} --snploc {snploc} --covariates {covariates} --keys {keys} --output-prefix {output_prefix} ',
        worker_boot_disk_size=200,
    )
    job.memory('standard')
    batch.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
