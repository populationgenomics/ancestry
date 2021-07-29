"""Entry point for the analysis runner."""

import os
import sys
import hailtop.batch as hb
from analysis_runner import dataproc


POP = sys.argv[1] if len(sys.argv) > 1 else 'nfe'

service_backend = hb.ServiceBackend(
    billing_project=os.getenv('HAIL_BILLING_PROJECT'), bucket=os.getenv('HAIL_BUCKET')
)

batch = hb.Batch(name=f'{POP} pca', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    f'hgdp_1kg_tob_wgs_pop_pca_densified.py --pop {POP}',
    max_age='4h',
    num_secondary_workers=20,
    packages=['click'],
    init=['gs://cpg-reference/hail_dataproc/install_common.sh'],
    job_name=f'{POP}-pca-new-variants',
    worker_boot_disk_size=200,
)

batch.run()
