"""Entry point for the analysis runner."""

import os
import hailtop.batch as hb
from analysis_runner import dataproc

service_backend = hb.ServiceBackend(
    billing_project=os.getenv('HAIL_BILLING_PROJECT'), bucket=os.getenv('HAIL_BUCKET')
)

batch = hb.Batch(name='plot-loadings-nfe-no-outliers', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    f'hgdp_1kg_tob_wgs_plot_loadings_nfe_no_outliers.py',
    max_age='4h',
    num_secondary_workers=20,
    packages=['selenium'],
    init=['gs://cpg-common-main/hail_dataproc/install_common.sh'],
    job_name=f'plot-loadings-nfe-no-outliers',
)

batch.run()
