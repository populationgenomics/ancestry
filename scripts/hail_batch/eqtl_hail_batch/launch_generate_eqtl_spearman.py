"""Launch analysis runner for all cell types and chromosomes"""

import re
import os
from analysis_runner.cli_analysisrunner import run_analysis_runner
from google.cloud import storage

REFERENCE = 'GRCh38'
BUCKET_NAME = 'cpg-tob-wgs-test'
INPUT_PATH = f'gs://{BUCKET_NAME}/full_data/{REFERENCE.lower()}/analysis_output/kccg'


# get cell types
cs_client = storage.Client()
bucket = cs_client.get_bucket(BUCKET_NAME)
bucket_path = INPUT_PATH.split(f'gs://{BUCKET_NAME}/')[-1]
blobs = [
    f'gs://{BUCKET_NAME}/{b.name}'
    for b in bucket.list_blobs(prefix=bucket_path, delimiter='/')
    if b.name.endswith('.csv')
]
cell_type_files = [
    str(path) for path in blobs if re.search(r'log_residuals', str(path))
]
cell_types = set(
    map(lambda filename: os.path.basename(filename).split('_')[0], cell_type_files)
)
n_chromosomes = list(range(1, 23))

for cell_type in cell_types:
    expression = f'gs://cpg-tob-wgs-test/kat/input/{cell_type}_expression.tsv'
    covariates = f'gs://cpg-tob-wgs-test/kat/input/{cell_type}_peer_factors.tsv'
    for idx in n_chromosomes:
        genotype = f'gs://cpg-tob-wgs-test/kat/input/genotype_chr{idx}.tsv'
        geneloc = f'gs://cpg-tob-wgs-test/kat/input/geneloc_chr{idx}.tsv'
        snploc = f'gs://cpg-tob-wgs-test/kat/input/snpsloc_chr{idx}.tsv'
        # if DRY_RUN:
        #     check_files_exists([expression, covariates, genotype, geneloc, snploc])
        # else:
        run_analysis_runner(
            script=[
                '--dataset',
                'tob-wgs',
                '--access-level',
                'test',
                '--output-dir',
                'kat/v0',
                '--description',
                'eqtl batch job',
                '--expression',
                expression,
                '--covariates',
                covariates,
                '--genotype',
                genotype,
                '--geneloc',
                geneloc,
                '--snploc',
                snploc,
                '--output-prefix',
                f'gs://cpg-tob-wgs-test/kat/{cell_type}/chr{idx}_',
            ]
        )
