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

batch = hb.Batch(name='test variant selection', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    f'hgdp_1kg_tob_wgs_variant_selection_exploration_test.py --output={OUTPUT}',
    max_age='5h',
    packages=['click', 'selenium'],
    init=['gs://cpg-reference/hail_dataproc/install_phantomjs.sh'],
    job_name='test-variant-selection',
)

batch.run()