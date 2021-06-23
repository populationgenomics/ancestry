"""Run nfe_loadings_qc.py using the analysis runner."""

import os
import hail as hl
import hailtop.batch as hb
from analysis_runner import dataproc

OUTPUT = os.getenv('OUTPUT')
assert OUTPUT

hl.init(default_reference='GRCh38')

service_backend = hb.ServiceBackend(
    billing_project=os.getenv('HAIL_BILLING_PROJECT'), bucket=os.getenv('HAIL_BUCKET')
)

batch = hb.Batch(name='loadings qc', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    f'nfe_loadings_qc.py --output={OUTPUT}',
    max_age='1h' 'nfe_loadings_qc.py',
    packages=['click'],
    job_name='loadings_qc',
)

batch.run()
