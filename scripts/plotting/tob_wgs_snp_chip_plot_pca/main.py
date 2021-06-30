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

batch = hb.Batch(name='pca_combined_tob_snp_chip', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    f'plot_pca_tob_wgs_snp_chip_datasets.py --output={OUTPUT}',
    max_age='1h',
    packages=['click', 'selenium'],
    init=['gs://cpg-reference/hail_dataproc/install_phantomjs.sh'],
    job_name=f'pca_combined_tob_snp_chip',
)

batch.run()
