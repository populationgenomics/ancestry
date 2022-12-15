"""Entry point for the analysis runner."""

import os
import hailtop.batch as hb
from analysis_runner import dataproc


service_backend = hb.ServiceBackend(
    billing_project=os.getenv('HAIL_BILLING_PROJECT'), bucket=os.getenv('HAIL_BUCKET')
)

batch = hb.Batch(name=f'increase_partitions', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    'increase_snp_chip_partitions.py',
    max_age='2h',
    num_workers=20,
    init=['gs://cpg-common-main/references/hail_dataproc/install_common.sh'],
    job_name=f'increase_partitions',
)

batch.run()
