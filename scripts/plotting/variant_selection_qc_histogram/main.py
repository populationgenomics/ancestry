"""Entry point for the analysis runner."""

import os
import hailtop.batch as hb
from analysis_runner import dataproc

service_backend = hb.ServiceBackend(
    billing_project=os.getenv('HAIL_BILLING_PROJECT'), bucket=os.getenv('HAIL_BUCKET')
)

batch = hb.Batch(name='variant selection', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    'variant_selection_qc_histogram.py',
    max_age='3h',
    packages=['click', 'selenium'],
    init=['gs://cpg-common-main/hail_dataproc/install_common.sh'],
    job_name='variant-selection-histogram',
)

batch.run()
