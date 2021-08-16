"""Entry point for the analysis runner."""

import os
import hailtop.batch as hb
from analysis_runner import dataproc

service_backend = hb.ServiceBackend(
    billing_project=os.getenv('HAIL_BILLING_PROJECT'), bucket=os.getenv('HAIL_BUCKET')
)

batch = hb.Batch(name='nfe-pc-relate', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    f'hgdp_1kg_tob_wgs_nfe_pc_relate.py',
    max_age='4h',
    num_secondary_workers=20,
    init=['gs://cpg-reference/hail_dataproc/install_common.sh'],
    job_name=f'nfe-pc-relate',
    worker_boot_disk_size=200,
)

batch.run()