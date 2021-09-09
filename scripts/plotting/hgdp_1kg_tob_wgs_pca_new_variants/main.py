"""Entry point for the analysis runner."""

import os
import hailtop.batch as hb
from analysis_runner import dataproc

service_backend = hb.ServiceBackend(
    billing_project=os.getenv('HAIL_BILLING_PROJECT'), bucket=os.getenv('HAIL_BUCKET')
)

batch = hb.Batch(name='subset-gnomad-hgdp-1kg', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    'subset_hgdp.py',
    max_age='12h',
    packages=['selenium'],
    init=['gs://cpg-reference/hail_dataproc/install_common.sh'],
    job_name='subset-gnomad-hgdp-1kg',
)

batch.run()
