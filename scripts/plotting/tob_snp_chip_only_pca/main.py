"""Entry point for the analysis runner."""

import os
import hailtop.batch as hb
from analysis_runner import dataproc

service_backend = hb.ServiceBackend(
    billing_project=os.getenv('HAIL_BILLING_PROJECT'), bucket=os.getenv('HAIL_BUCKET')
)

batch = hb.Batch(name='plot_snp_chip_pca', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    'plot_tob_snp_chip_pca_only.py',
    max_age='1h',
    packages=['selenium'],
    init=['gs://cpg-common-main/hail_dataproc/install_common.sh'],
    job_name=f'plot_snp_chip_pca',
)

batch.run()
