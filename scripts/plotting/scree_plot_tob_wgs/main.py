"""Entry point for the analysis runner."""

import os
import hailtop.batch as hb
from analysis_runner import dataproc

service_backend = hb.ServiceBackend(
    billing_project=os.getenv('HAIL_BILLING_PROJECT'), bucket=os.getenv('HAIL_BUCKET')
)

batch = hb.Batch(name='tob-wgs-scree', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    'tob_wgs_scree_plot.py',
    max_age='3h',
    packages=['selenium'],
    init=['gs://cpg-common-main/references/hail_dataproc/install_common.sh'],
    job_name='tob-wgs-scree',
)

batch.run()
