"""Entry point for the analysis runner."""

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

batch = hb.Batch(name=f'densify_tobwgs', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    f'hgdp_1kg_tob_wgs_densify.py --output={OUTPUT}',
    max_age='12h',
    num_secondary_workers=20,
    packages=['click'],
    job_name=f'densify_tobwgs',
)

batch.run()
