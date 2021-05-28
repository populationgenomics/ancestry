"""Run hgdp_1kg_tob_wgs_variant_selection.py using the analysis runner."""

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

batch = hb.Batch(name='variant selection', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    f'hgdp_1kg_tob_wgs_variant_selection.py --output={OUTPUT}',
    max_age='5h',
    num_secondary_workers=100,
    packages=['click'],
    job_name='variant-selection',
    worker_boot_disk_size=200,
)

batch.run()
