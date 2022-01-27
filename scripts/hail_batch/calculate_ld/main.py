"""Entry point for the analysis runner."""

import os
import hailtop.batch as hb
from analysis_runner import dataproc

service_backend = hb.ServiceBackend(
    billing_project=os.getenv('HAIL_BILLING_PROJECT'), bucket=os.getenv('HAIL_BUCKET')
)

batch = hb.Batch(name='calculate-ld', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    f'ld_exploration.py',
    max_age='12h',
    num_workers=20,
    init=['gs://cpg-reference/hail_dataproc/install_common.sh'],
    job_name=f'calculate_ld',
    worker_boot_disk_size=200,
)

batch.run()
