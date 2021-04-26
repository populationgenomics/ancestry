"""Run gnomad_loadings_90k_liftover.py using the analysis runner."""

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

batch = hb.Batch(name='gnomad loadings liftover', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    f'gnomad_loadings_90k_liftover.py --output={OUTPUT}',
    max_age='1h',
    packages=['click'],
    job_name='gnomad-loadings-liftover',
)

batch.run()
