"""Entry point for the analysis runner."""

import os
import hailtop.batch as hb
from analysis_runner import dataproc


service_backend = hb.ServiceBackend(
    billing_project=os.getenv('HAIL_BILLING_PROJECT'), bucket=os.getenv('HAIL_BUCKET')
)

batch = hb.Batch(name=f'project_wgs_samples', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    'project_wgs_samples_onto_snp_chip.py',
    max_age='12h',
    num_secondary_workers=20,
    packages=['selenium'],
    init=['gs://cpg-reference/hail_dataproc/install_common.sh'],
    job_name=f'project_wgs_samples',
)

batch.run()
