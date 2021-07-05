"""Entry point for the analysis runner."""

import os
import sys
import hail as hl
import hailtop.batch as hb
from analysis_runner import dataproc

OUTPUT = os.getenv('OUTPUT')
assert OUTPUT

hl.init(default_reference='GRCh38')

POP = sys.argv[1] if len(sys.argv) > 1 else 'nfe'

service_backend = hb.ServiceBackend(
    billing_project=os.getenv('HAIL_BILLING_PROJECT'), bucket=os.getenv('HAIL_BUCKET')
)

batch = hb.Batch(name=f'{POP} pc-project-test', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    f'pc_project_test.py --output={OUTPUT} --pop {POP}',
    max_age='5h',
    packages=['click', 'selenium'],
    job_name=f'{POP}-pc-project-test',
)

batch.run()
