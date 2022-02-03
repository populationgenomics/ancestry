"""Launch analysis runner for all cell types and chromosomes"""

import os
import click
from analysis_runner.cli_analysisrunner import run_analysis_runner
from google.cloud import storage


@click.command()
@click.option(
    '--cell_types',
    default=None,
    help='List of cell types to test',
)
@click.option(
    '--chromosomes_to_test',
    required=True,
    help='List of chromosome numbers to run eQTL analysis on.',
)
@click.option(
    '--input-path',
    required=True,
    help='A path prefix of where input files are located, eg: gs://MyBucket/folder/',
)
@click.option(
    '--output-dir',
    required=True,
    help='A path of where to output files, eg: gs://MyBucket/output-folder/',
)
def submit_eqtl_jobs(cell_types, chromosomes_to_test, input_path, output_dir):
    """Run association script for all chromosomes and cell types"""

    if cell_types is None:
        cs_client = storage.Client()
        bucket_name = input_path.split('gs://')[1].split('/')[0]
        bucket = cs_client.get_bucket(bucket_name)
        bucket_path = input_path.split(f'gs://{bucket_name}/')[-1]
        blobs = [
            f'gs://{bucket_name}/{b.name}'
            for b in bucket.list_blobs(prefix=bucket_path, delimiter='/')
            if b.name.endswith('.tsv')
        ]
        # get celltype name, which is the first word in the filename
        cell_types = set(os.path.basename(filename).split('_')[0] for filename in blobs)

    for cell_type in cell_types:
        expression = f'gs://cpg-tob-wgs-test/kat/input/{cell_type}_expression.tsv'
        covariates = f'gs://cpg-tob-wgs-test/kat/input/{cell_type}_peer_factors.tsv'
        for idx in chromosomes_to_test:
            genotype = f'gs://cpg-tob-wgs-test/kat/input/genotype_chr{idx}.tsv'
            geneloc = f'gs://cpg-tob-wgs-test/kat/input/geneloc_chr{idx}.tsv'
            snploc = f'gs://cpg-tob-wgs-test/kat/input/snpsloc_chr{idx}.tsv'
            # if DRY_RUN: check that all the files exist
            #     check_files_exists([expression, covariates, 
            # genotype, geneloc, snploc])
            # else:
            run_analysis_runner(
                script=[
                    '--dataset',
                    'tob-wgs',
                    '--access-level',
                    'test',
                    '--output-dir',
                    output_dir,
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
                    f'gs://cpg-tob-wgs-test/kat/{output_dir}',
                ]
            )


if __name__ == '__main__':
    submit_eqtl_jobs()
