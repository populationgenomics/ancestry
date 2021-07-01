"""Entry point for the analysis runner."""

import os
import hailtop.batch as hb
from analysis_runner import dataproc


service_backend = hb.ServiceBackend(
    billing_project=os.getenv('HAIL_BILLING_PROJECT'), bucket=os.getenv('HAIL_BUCKET')
)

batch = hb.Batch(name=f'snp_chip_variants_pca', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    'snp_chip_generate_pca.py',
    max_age='12h',
    num_secondary_workers=20,
    packages=['click'],
    init=['gs://cpg-reference/hail_dataproc/install_common.sh'],
    job_name=f'snp_chip_variants_pca',
)

batch.run()
