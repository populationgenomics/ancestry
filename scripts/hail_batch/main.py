"""Run hgdp_1kg_pca_loadings_generate.py using the analysis runner."""

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

batch = hb.Batch(name='generate PCA loadings', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    f'HGDP_1KG_PCALoadingsGenerate.py --output={OUTPUT}',
    max_age='15h',
    num_secondary_workers='100',
    packages=['click'],
    job_name='PCA-loadings',
)

batch.run()