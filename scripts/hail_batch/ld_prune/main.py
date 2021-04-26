"""Run hgdp_1kg_ld_prune.py using the analysis runner."""

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

batch = hb.Batch(name='ld pruning', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    f'hgdp_1kg_ld_prune.py --output={OUTPUT}',
    max_age='5h',
    num_secondary_workers=300,
    packages=['click', 'gnomad'],
    job_name='ld-prune',
)

batch.run()
