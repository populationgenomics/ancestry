"""Entry point for the analysis runner."""

import os
import hailtop.batch as hb
from analysis_runner import dataproc


service_backend = hb.ServiceBackend(
    billing_project=os.getenv('HAIL_BILLING_PROJECT'), bucket=os.getenv('HAIL_BUCKET')
)

batch = hb.Batch(name=f'densify_tobwgs_new_variants', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    'hgdp_1kg_tob_wgs_densify_new_variants.py',
    max_age='12h',
    num_secondary_workers=20,
    packages=['click'],
    init=['gs://cpg-reference/hail_dataproc/install_common.sh'],
    job_name=f'densify_tobwgs_pca',
)

batch.run()
