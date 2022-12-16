"""Entry point for the analysis runner."""

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

batch = hb.Batch(name='densified pca nfe', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    f'hgdp_1kg_tob_wgs_plot_pca_nfe.py --output={OUTPUT}',
    max_age='1h',
    packages=['click', 'selenium'],
    init=['gs://cpg-common-main/hail_dataproc/install_phantomjs.sh'],
    job_name=f'densified pca nfe',
)

batch.run()
