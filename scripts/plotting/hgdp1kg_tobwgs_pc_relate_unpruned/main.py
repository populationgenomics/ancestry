"""Entry point for the analysis runner."""

import os
import hailtop.batch as hb
from analysis_runner import dataproc

service_backend = hb.ServiceBackend(
    billing_project=os.getenv('HAIL_BILLING_PROJECT'), bucket=os.getenv('HAIL_BUCKET')
)

batch = hb.Batch(name='related_samples-plot', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    f'hgdp_1kg_tob_wgs_related_samples.py',
    max_age='2h',
    init=['gs://cpg-reference/hail_dataproc/install_common.sh'],
    job_name=f'related_samples-plot',
)

batch.run()
