#!/usr/bin/env python3

"""
Use VEP using a dataproc cluster.
Taken from Matt Welland's script, run_vep_help.py

"""


import os
import hailtop.batch as hb
from analysis_runner import dataproc
import click


@click.command()
@click.option('--script', 'script', help='path to VEP main script')
def main(script: str):
    """
    runs a script inside dataproc to execute VEP
    :param script: str, the path to the VEP main script
    """

    service_backend = hb.ServiceBackend(
        billing_project=os.getenv('HAIL_BILLING_PROJECT'),
        bucket=os.getenv('HAIL_BUCKET'),
    )

    # create a hail batch
    batch = hb.Batch(name='run_vep_in_dataproc_cluster', backend=service_backend)

    job = dataproc.hail_dataproc_job(
        batch=batch,
        script=script,
        max_age='12h',
        job_name='run_vep',
        num_secondary_workers=20,
        cluster_name='run vep',
        vep='GRCh38',
    )
    job.cpu(2)
    job.memory('standard')
    job.storage('20G')

    batch.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
