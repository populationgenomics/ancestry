"""Perform conditional analysis on SNPs and expression residuals"""

import os
import hail as hl
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


# Run Spearman rank in parallel by sending genes in a batches
def run_computation_in_scatter(idx, inputs=None):  # pylint: disable=too-many-locals
    """Run genes in scatter"""
    # Input filenames
    if inputs is None:
        residual_df = pd.read_csv(
            f'gs://{INPUT_BUCKET}/kat/input/Plasma_chr22_log_residuals.tsv', sep='\t'
        )
        significant_snps = pd.read_csv(
            f'gs://{INPUT_BUCKET}/kat/test.csv', sep=' ', skipinitialspace=True
        )
    else:
        # if you specify the inputs, they're in this order (so deconstruct them)
        residual_df, significant_snps = inputs

    genotype_df = pd.read_csv(
        f'gs://{INPUT_BUCKET}/kat/input/genotype_chr22.tsv', sep='\t'
    )

    # Identify the top eSNP for each eGene and assign remaining to df
    esnp1 = (
        significant_snps.sort_values(['geneid', 'p.value'], ascending=True)
        .groupby('geneid')
        .first()
        .reset_index()
    )
    print('Printing esnp1')
    print(esnp1)
    print('Printing significant snps')
    print(significant_snps)
    esnps_to_test = (
        significant_snps.sort_values(['geneid', 'p.value'], ascending=True)
        .groupby('geneid')
        .apply(lambda group: group.iloc[1:, 1:])
        # .reset_index()
    )
    print(esnps_to_test)
    esnps_to_test = esnps_to_test.reset_index()

    # Subset residuals for the genes to be tested
    sample_ids = residual_df.loc[:, ['sampleid']]
    gene_ids = esnp1['geneid'][esnp1['geneid'].isin(residual_df.columns)]
    residual_df = residual_df.loc[:, residual_df.columns.isin(gene_ids)]
    residual_df['sampleid'] = sample_ids

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

    # Spearman's rank correlation
    def spearman_correlation(df):
        """get Spearman rank correlation"""
        gene = df.geneid
        snp = df.snpid
        res_val = adjusted_residual_mat[['sampleid', gene]]
        genotype_val = genotype_df[['sampleid', snp]]
        test_df = res_val.merge(genotype_val, on='sampleid', how='left')
        test_df.columns = ['sampleid', 'residual', 'SNP']
        coef, p = spearmanr(test_df['SNP'], test_df['residual'])
        return (gene, snp, coef, p)

    esnps_to_test = esnps_to_test[
        esnps_to_test.geneid.isin(adjusted_residual_mat.columns)
    ]
    gene_snp_test_df = esnps_to_test[['snpid', 'geneid']]
    gene_snp_test_df = gene_snp_test_df[
        gene_snp_test_df['geneid'] == gene_ids.iloc[idx]
    ]
    adjusted_spearman_df = pd.DataFrame(
        list(gene_snp_test_df.apply(spearman_correlation, axis=1))
    )
    adjusted_spearman_df.columns = ['geneid', 'snpid', 'coef', 'p.value']
    # add in global position and round
    locus = adjusted_spearman_df.snpid.str.split('_', expand=True)[0]
    chromosome = 'chr' + locus.str.split(':', expand=True)[0]
    position = locus.str.split(':', expand=True)[1]
    (
        adjusted_spearman_df['locus'],
        adjusted_spearman_df['chromosome'],
        adjusted_spearman_df['position'],
    ) = [
        locus,
        chromosome,
        position,
    ]
    adjusted_spearman_df['round'] = 2
    # convert to hail table. Can't call `hl.from_pandas(spearman_df)` directly
    # because it doesnt' work with the spark local backend
    adjusted_spearman_df.to_csv(f'adjusted_spearman_df.csv')
    hl.init(default_reference='GRCh38')
    t = hl.import_table(
        'adjusted_spearman_df.csv',
        delimiter=',',
        types={'position': hl.tint32, 'coef': hl.tfloat64, 'p.value': hl.tfloat64},
    )
    t = t.annotate(global_position=hl.locus(t.chromosome, t.position).global_position())
    # get alleles
    # mt.rows()[c.liftover].alleles
    # turn back into pandas df. Can't call `spearman_df = t.to_pandas()` directly
    # because it doesn't work with the spark local backend
    t.export('adjusted_spearman_df_annotated.tsv')
    adjusted_spearman_df = pd.read_csv('adjusted_spearman_df_annotated.tsv', sep='\t')

    # set variables for next iteration of loop
    residual_df = adjusted_residual_mat
    significant_snps = adjusted_spearman_df

    # the order of this is important
    return [residual_df, significant_snps]


def function_that_merges_lists_of_dataframes(*df_list):
    """
    Merge list of list of dataframes
    df_list = [
        (residual_df_1, significant_snps_1),
        (residual_df_2, significant_snps_2),
        ...
    ]
    returns [merged_residual_df, merged_significant_snps]

    """
    merged_dfs = []
    for idx in range(len(df_list[0])):
        merged_df = pd.concat([i[idx] for i in df_list])
        pvalues = merged_df['p.value']
        fdr_values = pd.DataFrame(list(multi.fdrcorrection(pvalues))).iloc[1]
        merged_df = merged_df.assign(FDR=fdr_values)
        merged_df['FDR'] = merged_df.FDR.astype(float)
        merged_dfs.append(merged_df)
    return merged_dfs


backend = hb.ServiceBackend(billing_project='tob-wgs', bucket='cpg-tob-wgs-test')
b = hb.Batch(name='eQTL', backend=backend, default_python_image=DRIVER_IMAGE)

# significants_snps = []
# iteratively, do 5 times:
#       For each gene:
#           find residuals
#           find spearman rank correlation
#       merge(residuals)
#       merge significant_snps
#       ^ feed those two into next iteration

N_GENES = 5
# N_GENES = get_number_of_scatters()
# for i in range(get_number_of_scatters()):

result = None  # pylint: disable=invalid-name
for _ in range(5):
    residual_and_sig_snps_dfs = []
    for i in range(N_GENES):
        j = b.new_python_job(name=f'process_{i}')
        result: hb.resource.PythonResult = j.call(run_computation_in_scatter, i, result)
        residual_and_sig_snps_dfs.append(result)

    merge_job = b.new_python_job(name='merge_scatters')
    result = merge_job.call(
        function_that_merges_lists_of_dataframes, *residual_and_sig_snps_dfs
    )


def get_last_element_of_list_and_convert_to_text(element):
    """
    the result array is of format [residual_df, significant_snps]
    so get LAST element and convert to string
    """
    return element[-1].to_string()


get_sig_snps_job = b.new_python_job('get_significant_snps')
sig_snps_as_string = get_sig_snps_job.call(
    get_last_element_of_list_and_convert_to_text, result
)

b.write_output(
    sig_snps_as_string.as_str(),
    f'gs://{OUTPUT_BUCKET}/kat/test_conditional_analysis.csv',
)
b.run()
