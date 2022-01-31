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
blobs = [
    f'gs://{BUCKET_NAME}/{b.name}'
    for b in bucket.list_blobs()
    if (INPUT_PATH.split(f'gs://{BUCKET_NAME}/')[-1] in b.name) and ('.tsv' in b.name)
]
cell_type_files = [
    str(path) for path in blobs if re.search(r'log_residuals', str(path))
]
cell_types = set(
    map(lambda filename: os.path.basename(filename).split('_')[0], cell_type_files)
)
n_chromosomes = list(range(1, 23))

for cell_type in cell_types:
    for idx in n_chromosomes:
        residuals = (
            f'gs://cpg-tob-wgs-test/kat/input/{cell_type}/chr{idx}_log_residuals.tsv'
        )
        genotype = f'gs://cpg-tob-wgs-test/kat/input/genotype_chr{idx}.tsv'
        significant_snps = (
            f'gs://cpg-tob-wgs-test/kat/input/{cell_type}/'
            f'chr{idx}_correlation_results.tsv'
        )
        # if DRY_RUN:
        #     check_files_exists([residuals, genotype, significant_snps])
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
                '--residuals',
                residuals,
                '--significant_snps',
                significant_snps,
                '--genotype',
                genotype,
                '--output-prefix',
                f'gs://cpg-tob-wgs-test/kat/{cell_type}/chr{idx}_'
                '--test_subset_genes',
                '5',
            ]
        )
