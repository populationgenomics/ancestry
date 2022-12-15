"""Entry point for the analysis runner."""

import os
import hailtop.batch as hb
from analysis_runner import dataproc

service_backend = hb.ServiceBackend(
    billing_project=os.getenv('HAIL_BILLING_PROJECT'), bucket=os.getenv('HAIL_BUCKET')
)

batch = hb.Batch(name=f'tobwgs_pca_new_variants', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    'hgdp_1kg_tob_wgs_densified_pca_new_variants.py',
    max_age='12h',
    num_secondary_workers=20,
    init=['gs://cpg-common-main/references/hail_dataproc/install_common.sh'],
    job_name=f'tobwgs_pca_new_variants',
    worker_boot_disk_size=200,
)

batch.run()
