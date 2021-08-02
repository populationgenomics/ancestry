"""Entry point for the analysis runner."""

import os
import hailtop.batch as hb
from analysis_runner import dataproc

service_backend = hb.ServiceBackend(
    billing_project=os.getenv('HAIL_BILLING_PROJECT'), bucket=os.getenv('HAIL_BUCKET')
)

batch = hb.Batch(name='no-outliers-plot-pca-nfe', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    f'hgdp_1kg_tob_wgs_plot_pca_nfe.py',
    max_age='1h',
    packages=['selenium'],
    init=['gs://cpg-reference/hail_dataproc/install_common.sh'],
    job_name='no-outliers-plot-pca-nfe',
)

batch.run()
