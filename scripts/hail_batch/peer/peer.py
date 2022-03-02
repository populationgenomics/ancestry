"""
Test run this with:

    cd <dir where peer.py file is>
    analysis-runner \
        --access-level test \
        --description 'Test run peer analysis' \
        --output-dir 'kat/2022-03-2_peer' \
        --dataset tob-wgs \
        python peer.py
"""


import os
import hailtop.batch as hb
import pandas as pd

DRIVER_IMAGE = 'australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:e6451763492b62ddfadc20c06b240234b20b6f2f-hail-0.2.73.devc6f6f09cec08'
PEER_DOCKER = 'australia-southeast1-docker.pkg.dev/cpg-common/images/peer:1.3.1'

SCORES_PATH = 'gs://cpg-tob-wgs-test/kat/pca/nfe_feb22/v0/scores.json'
COVARIATES_PATH = 'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/covariates_files/covariates.tsv'
SAMPLE_ID_KEYS_PATH = (
    'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/OneK1K_CPG_IDs.tsv'
)


def get_covariates(scores_path, covariates_path, sample_id_keys_path) -> str:
    """
    Get covariate data by merging PCA scores with age and sex info.
    Only needs to be run once. This returns a TSV (as a string)
    """
    scores_df = pd.read_json(scores_path)
    sampleid = scores_df.s
    # only get the first 4 PCs, as indicated by scree plot
    scores = pd.DataFrame(scores_df['scores'].to_list()).iloc[:, 0:4]
    scores.columns = ['PC1', 'PC2', 'PC3', 'PC4']
    scores.insert(loc=0, column='sampleid', value=sampleid)

    # get age and sex from Seyhan's covariate file
    covariates = pd.read_csv(covariates_path, sep='\t')
    # load in keys file to convert sampleids in covariates file to CPG IDs
    sampleid_keys = pd.read_csv(sample_id_keys_path, sep='\t')
    covariates = pd.merge(
        covariates, sampleid_keys, how='left', left_on='sampleid', right_on='OneK1K_ID'
    ).drop(['pc1', 'pc2', 'pc3', 'pc4'], axis=1)
    # merge together covariates df and scores
    covariates = pd.merge(
        covariates, scores, how='left', left_on='InternalID', right_on='sampleid'
    ).drop(['OneK1K_ID', 'InternalID', 'ExternalID', 'sampleid_y'], axis=1)

    return covariates.to_csv()


def run_peer_job(b: hb.Batch, expression_file, covariates_file):
    """
    Run peer analysis, except because peer is in Python2, we can't take
    advantage of the python_job concept in Batch. Instead, we'll rely on
    the two inputs `expression_file` and `covariates_file` to be TSV files.
    """

    j = b.new_job('peer')

    j.image(PEER_DOCKER)

    j.command(f'cat {expression_file}')
    j.command(f'cat {covariates_file}')

    # write python script to container
    j.command(
        """
cat <<EOT >> run_peer.py

import peer
import numpy as np
    
def run_peer(expression_file, covariates_file, factors_output_path):
    \"""
    Get covariate data for each cell type
    \"""

    print 'Loading data'

    # load in data
    expr = np.loadtxt(expression_file, delimiter=',')
    covs = np.loadtxt(covariates_file, delimiter=',')

    print 'Loaded data'

    # Set PEER paramaters as per the PEER website
    model = peer.PEER()
    model.setPhenoMean(expr)
    model.setCovariates(covs)
    model.setNk(10)
    model.update()

    # Calculate and save the PEER factors
    factors = model.getX()
    np.savetxt(factors_output_path, factors, delimiter=',')
    # Calculate and save the weights for each factor
    weights = model.getW()
    np.savetxt('weight.csv', weights, delimiter=',')
    # Calculate and save the precision values
    precision = model.getAlpha()
    np.savetxt('precision.csv', precision, delimiter=',')
    # Calculate and save the residuals
    residuals = model.getResiduals()
    np.savetxt('residuals.csv', residuals, delimiter=',')


if __name__ == "__main__":
    import sys
    run_peer(sys.argv[1], sys.argv[2], sys.argv[3])
    
EOT"""
    )

    j.command(
        f'python run_peer.py {expression_file} {covariates_file} {j.factors_output_path}'
    )
    j.command('ls -l')

    return j


# @click.command()
# @click.option('--expression-file')
def main(
    expression_file,
    scores_path=SCORES_PATH,
    covariates_path=COVARIATES_PATH,
    sample_id_keys_path=SAMPLE_ID_KEYS_PATH,
):
    """
    Run PEER calculation in hail batch
    """

    dataset = os.getenv('DATASET', 'tob-wgs')
    access_level = os.getenv('ACCESS_LEVEL', 'test')
    driver_image = os.getenv('DRIVER_IMAGE', DRIVER_IMAGE)

    backend = hb.ServiceBackend(
        billing_project=dataset, bucket=f'cpg-{dataset}-{access_level}'
    )

    batch = hb.Batch(name='PEER', backend=backend, default_python_image=driver_image)
    expression_f = batch.read_input(expression_file)

    load_covariates = batch.new_python_job('load covariates')

    covariates_tsv = load_covariates.call(
        get_covariates, scores_path, covariates_path, sample_id_keys_path
    ).as_str()

    peer_job = run_peer_job(batch, expression_f, covariates_tsv)

    second_job = batch.new_python_job('downstream tasks')
    second_job.call(print, peer_job.factors_output_path)

    # batch.run(dry_run=True)
    batch.run(wait=False)


if __name__ == '__main__':
    main(
        expression_file='gs://cpg-tob-wgs-test/kat/expression.csv',
    )  # pylint: disable=no-value-for-parameter
