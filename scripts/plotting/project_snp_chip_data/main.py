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

batch = hb.Batch(name=f'{POP} project snp-chip', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    f'project_snp_chip_data.py --output={OUTPUT} --pop {POP}',
    max_age='5h',
    packages=['click', 'selenium'],
    init=['gs://cpg-common-main/references/hail_dataproc/install_phantomjs.sh'],
    job_name=f'{POP}-project-snp-chip',
)

batch.run()
