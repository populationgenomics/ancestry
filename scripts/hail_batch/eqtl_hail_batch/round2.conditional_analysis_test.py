"""Perform conditional analysis on SNPs and expression residuals"""

import os

# import hail as hl
import hailtop.batch as hb
import pandas as pd
import statsmodels.api as sm
import statsmodels.stats.multitest as multi
from patsy import dmatrices  # pylint: disable=no-name-in-module
from scipy.stats import spearmanr

INPUT_BUCKET = 'cpg-tob-wgs-test'
OUTPUT_BUCKET = 'cpg-tob-wgs-test'
DRIVER_IMAGE = os.getenv(
    'DRIVER_IMAGE',
    'australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d2a9c316d6d752edb27623542c8a062db4466842-hail-0.2.73.devc6f6f09cec08',  # noqa: E501; pylint: disable=line-too-long
)
DEFAULT_RESIDUALS_PATH = f'gs://{INPUT_BUCKET}/kat/input/Plasma_chr22_log_residuals.tsv'
DEFAULT_SNPS_PATH = f'gs://{INPUT_BUCKET}/kat/correlation_results.csv'

# def get_number_of_scatters():
#     """get index of total number of genes"""
#     # Input filenames
#     residual_df = pd.read_csv(
#         f'gs://{INPUT_BUCKET}/kat/input/Plasma_chr22_log_residuals.tsv', sep='\t'
#     )
#     significant_snps = pd.read_csv(
#         f'gs://{INPUT_BUCKET}/kat/test.csv', sep=' ', skipinitialspace=True
#     )

#     # Identify the top eSNP for each eGene and assign remaining to df
#     esnp1 = (
#         significant_snps.sort_values(['geneid', 'p.value'], ascending=True)
#         .groupby('geneid')
#         .first()
#         .reset_index()
#     )
#     gene_ids = esnp1['geneid'][esnp1['geneid'].isin(residual_df.columns)]

#     return len(gene_ids)


def get_genotype_df(significant_snps, sample_ids):
    """load genotype df and filter"""
    genotype_df = pd.read_csv(
        f'gs://{INPUT_BUCKET}/kat/input/genotype_chr22.tsv', sep='\t'
    )
    # Subset genotype file for the significant SNPs
    genotype_df = genotype_df.loc[
        genotype_df['sampleid'].isin(  # pylint: disable=unsubscriptable-object
            sample_ids.sampleid
        ),
        :,
    ]
    genotype_sampleid = genotype_df['sampleid']
    genotype_df = genotype_df.loc[:, genotype_df.columns.isin(significant_snps.snpid)]
    genotype_df = genotype_df.assign(sampleid=genotype_sampleid)

    return genotype_df


def calculate_residual_df(residual_df, significant_snps):
    """calculate residuals for gene list"""
    if residual_df is None:
        residual_df = pd.read_csv(DEFAULT_RESIDUALS_PATH, sep='\t')
    if significant_snps is None:
        significant_snps = pd.read_csv(
            DEFAULT_SNPS_PATH, sep=' ', skipinitialspace=True
        )

    # Identify the top eSNP for each eGene and assign remaining to df
    esnp1 = (
        significant_snps.sort_values(['geneid', 'FDR'], ascending=True)
        .groupby('geneid')
        .first()
        .reset_index()
    )

    # Subset residuals for the genes to be tested
    sample_ids = residual_df.loc[:, ['sampleid']]
    gene_ids = esnp1['geneid'][esnp1['geneid'].isin(residual_df.columns)]
    residual_df = residual_df.loc[:, residual_df.columns.isin(gene_ids)]
    residual_df['sampleid'] = sample_ids
    genotype_df = get_genotype_df(significant_snps, sample_ids)

    # Find residuals after adjustment of lead SNP
    def calculate_adjusted_residuals(gene_id):
        gene = gene_id
        # select gene to regress
        exprs_val = residual_df[['sampleid', gene]]
        # select SNP to add
        snp = list(esnp1.snpid[esnp1.geneid == gene].values)
        snp.append('sampleid')
        snp_genotype = genotype_df[snp]

        # Create a test df by adding covariates
        test_df = exprs_val.merge(snp_genotype, on='sampleid', how='left')
        test_df.columns = ['sampleid', 'expression', 'genotype']

        y, x = dmatrices('expression ~ genotype', test_df)
        model = sm.OLS(y, x)
        residuals = list(model.fit().resid)
        return residuals

    adjusted_residual_mat = pd.DataFrame(
        list(map(calculate_adjusted_residuals, gene_ids))
    ).T
    adjusted_residual_mat.columns = gene_ids
    adjusted_residual_mat.insert(loc=0, column='sampleid', value=sample_ids.sampleid)
    # adjusted_residual_mat.to_csv(
    #     f'gs://{OUTPUT_BUCKET}/kat/eSNP1_adjusted_residuals.tsv'
    # )

    return adjusted_residual_mat


# Run Spearman rank in parallel by sending genes in a batches
def run_computation_in_scatter(
    iteration,  # pylint: disable=redefined-outer-name
    idx,
    residual_df,
    significant_snps=None,
):
    """Run genes in scatter"""
    # Input filenames
    if significant_snps is None:
        significant_snps = pd.read_csv(
            DEFAULT_SNPS_PATH, sep=' ', skipinitialspace=True
        )

    print(f'iteration = {iteration}')
    print(f'idx = {idx}')
    # make sure 'geneid' is the first column
    # otherwise, error thrown when using reset_index
    cols = list(significant_snps)
    cols.insert(0, cols.pop(cols.index('geneid')))
    significant_snps = significant_snps.loc[:, cols]
    esnps_to_test = (
        significant_snps.sort_values(['geneid', 'FDR'], ascending=True)
        .groupby('geneid')
        .apply(lambda group: group.iloc[1:, 1:])
        .reset_index(drop=True)
    )

    sample_ids = residual_df.loc[:, ['sampleid']]
    genotype_df = get_genotype_df(significant_snps, sample_ids)

    # Spearman's rank correlation
    def spearman_correlation(df):
        """get Spearman rank correlation"""
        gene = df.geneid
        snp = df.snpid
        res_val = residual_df[['sampleid', gene]]
        genotype_val = genotype_df[['sampleid', snp]]
        test_df = res_val.merge(genotype_val, on='sampleid', how='left')
        test_df.columns = ['sampleid', 'residual', 'SNP']
        coef, p = spearmanr(test_df['SNP'], test_df['residual'])
        return (gene, snp, coef, p)

    esnp1 = (
        significant_snps.sort_values(['geneid', 'FDR'], ascending=True)
        .groupby('geneid')
        .first()
        .reset_index()
    )
    gene_ids = esnp1['geneid'][esnp1['geneid'].isin(residual_df.columns)]

    esnps_to_test = esnps_to_test[esnps_to_test.geneid.isin(residual_df.columns)]
    gene_snp_test_df = esnps_to_test[['snpid', 'geneid']]
    gene_snp_test_df = gene_snp_test_df[
        gene_snp_test_df['geneid'] == gene_ids.iloc[idx]
    ]
    adjusted_spearman_df = pd.DataFrame(
        list(gene_snp_test_df.apply(spearman_correlation, axis=1))
    )
    adjusted_spearman_df.columns = ['geneid', 'snpid', 'coef', 'p.value']
    # add in global position and round
    # locus = adjusted_spearman_df.snpid.str.split('_', expand=True)[0]
    # chromosome = 'chr' + locus.str.split(':', expand=True)[0]
    # position = locus.str.split(':', expand=True)[1]
    # (
    #     adjusted_spearman_df['locus'],
    #     adjusted_spearman_df['chromosome'],
    #     adjusted_spearman_df['position'],
    # ) = [
    #     locus,
    #     chromosome,
    #     position,
    # ]
    adjusted_spearman_df['round'] = iteration
    # convert to hail table. Can't call `hl.from_pandas(spearman_df)` directly
    # because it doesnt' work with the spark local backend
    # adjusted_spearman_df.to_csv(f'adjusted_spearman_df.csv')
    # hl.init(default_reference='GRCh38')
    # t = hl.import_table(
    #     'adjusted_spearman_df.csv',
    #     delimiter=',',
    #     types={'position': hl.tint32, 'coef': hl.tfloat64, 'p.value': hl.tfloat64},
    # )
    # t = t.annotate(global_position=hl.locus(t.chromosome, t.position).global_position()) # noqa: E501; pylint: disable=line-too-long
    # get alleles
    # mt.rows()[c.liftover].alleles
    # turn back into pandas df. Can't call `spearman_df = t.to_pandas()` directly
    # because it doesn't work with the spark local backend
    # t.export('adjusted_spearman_df_annotated.tsv')
    # adjusted_spearman_df = pd.read_csv('adjusted_spearman_df_annotated.tsv', sep='\t')

    # set variables for next iteration of loop
    significant_snps = adjusted_spearman_df

    return significant_snps


def merge_significant_snps_dfs(*df_list):
    """
    Merge list of list of sig_snps dataframes
    """

    # merged sig_snps
    merged_sig_snps = pd.concat(df_list)
    pvalues = merged_sig_snps['p.value']
    fdr_values = pd.DataFrame(list(multi.fdrcorrection(pvalues))).iloc[1]
    merged_sig_snps = merged_sig_snps.assign(FDR=fdr_values)
    merged_sig_snps['FDR'] = merged_sig_snps.FDR.astype(float)
    merged_sig_snps.append(merged_sig_snps)

    return merged_sig_snps


def convert_dataframe_to_text(df):
    """
    convert to string
    """
    return df.to_string()


backend = hb.ServiceBackend(billing_project='tob-wgs', bucket='cpg-tob-wgs-test')
b = hb.Batch(name='eQTL', backend=backend, default_python_image=DRIVER_IMAGE)

N_GENES = 5
# N_GENES = get_number_of_scatters()
# for i in range(get_number_of_scatters()):

previous_sig_snps_result = None  # pylint: disable=invalid-name
previous_residual_result = None  # pylint: disable=invalid-name
for iteration in range(5):

    calc_resid_df_job = b.new_python_job(f'calculate-resid-df-iter-{iteration}')
    previous_residual_result = calc_resid_df_job.call(
        calculate_residual_df, previous_residual_result, previous_sig_snps_result
    )

    sig_snps_dfs = []
    for i in range(N_GENES):
        j = b.new_python_job(name=f'process_iter_{iteration}_job_{i}')
        gene_result: hb.resource.PythonResult = j.call(
            run_computation_in_scatter,
            iteration,
            i,
            previous_residual_result,
            previous_sig_snps_result,
        )
        sig_snps_dfs.append(gene_result)

    merge_job = b.new_python_job(name='merge_scatters')
    previous_sig_snps_result = merge_job.call(merge_significant_snps_dfs, *sig_snps_dfs)

    # convert sig snps to string for output
    sig_snps_as_string = merge_job.call(
        convert_dataframe_to_text, previous_sig_snps_result
    )
    # output sig snps for each iteration
    b.write_output(
        sig_snps_as_string.as_str(),
        f'gs://{OUTPUT_BUCKET}/kat/round{iteration+1}_significant_correlation_results.csv',  # noqa: E501; pylint: disable=line-too-long
    )

b.run(wait=False)
