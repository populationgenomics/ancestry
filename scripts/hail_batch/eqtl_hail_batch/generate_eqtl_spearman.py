"""Run Spearman rank correlation on SNPs and expression residuals"""

import os

# import hail as hl
import hailtop.batch as hb
import pandas as pd
import numpy as np
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


def get_number_of_scatters():
    """get index of total number of genes"""
    expression_df = pd.read_csv(
        f'gs://{INPUT_BUCKET}/kat/input/Plasma_expression.tsv', sep='\t'
    )
    # Remove genes with 0 expression in all samples
    expression_df = expression_df.loc[:, (expression_df != 0).any(axis=0)]
    # Number of individuals with non-zero expression
    genes_not_equal_zero = expression_df.iloc[:, 1:].values != 0
    n_expr_over_zero = pd.DataFrame(genes_not_equal_zero.sum(axis=0))
    percent_expr_over_zero = (n_expr_over_zero / len(expression_df.index)) * 100
    percent_expr_over_zero.index = expression_df.columns[1:]

    # Filter genes with less than 10 percent individuals with non-zero expression
    atleast10percent = percent_expr_over_zero[(percent_expr_over_zero > 10)[0]]
    sample_ids = expression_df['sampleid']
    expression_df = expression_df[atleast10percent.index]
    expression_df.insert(loc=0, column='sampleid', value=sample_ids)
    gene_ids = list(expression_df.columns.values)[1:]
    geneloc_df = pd.read_csv(
        f'gs://{INPUT_BUCKET}/kat/input/geneloc_chr22.tsv', sep='\t'
    )
    geneloc_df = geneloc_df[geneloc_df.geneid.isin(gene_ids)]

    return len(geneloc_df.index)


# Run Spearman rank in parallel by sending genes in a batches
def run_computation_in_scatter(idx):  # pylint: disable=too-many-locals
    """Run genes in scatter"""

    expression_df = pd.read_csv(
        f'gs://{INPUT_BUCKET}/kat/input/Plasma_expression.tsv', sep='\t'
    )
    genotype_df = pd.read_csv(
        f'gs://{INPUT_BUCKET}/kat/input/genotype_chr22.tsv', sep='\t'
    )
    geneloc_df = pd.read_csv(
        f'gs://{INPUT_BUCKET}/kat/input/geneloc_chr22.tsv', sep='\t'
    )
    snploc_df = pd.read_csv(
        f'gs://{INPUT_BUCKET}/kat/input/snpsloc_chr22.tsv', sep='\t'
    )
    covariate_df = pd.read_csv(
        f'gs://{INPUT_BUCKET}/kat/input/Plasma_peer_factors.tsv', sep='\t'
    )

    # Remove genes with 0 expression in all samples
    expression_df = expression_df.loc[:, (expression_df != 0).any(axis=0)]
    # Number of individuals with non-zero expression
    genes_not_equal_zero = expression_df.iloc[:, 1:].values != 0
    n_expr_over_zero = pd.DataFrame(genes_not_equal_zero.sum(axis=0))
    percent_expr_over_zero = (n_expr_over_zero / len(expression_df.index)) * 100
    percent_expr_over_zero.index = expression_df.columns[1:]

    # Filter genes with less than 10 percent individuals with non-zero expression
    atleast10percent = percent_expr_over_zero[(percent_expr_over_zero > 10)[0]]
    sample_ids = expression_df['sampleid']
    expression_df = expression_df[atleast10percent.index]
    expression_df.insert(loc=0, column='sampleid', value=sample_ids)

    # Prepare variables used to calculate Spearman's correlation
    gene_ids = list(expression_df.columns.values)[1:]
    sample_ids = expression_df.iloc[:, 0]
    genotype_df = genotype_df[genotype_df.sampleid.isin(expression_df.sampleid)]

    # Get 1Mb sliding window around each gene
    geneloc_df = geneloc_df[geneloc_df.geneid.isin(gene_ids)]
    geneloc_df = geneloc_df.assign(left=geneloc_df.start - 1000000)
    geneloc_df = geneloc_df.assign(right=geneloc_df.end + 1000000)

    to_log = expression_df.iloc[:, 1:].columns
    log_expression_df = expression_df[to_log].applymap(lambda x: np.log(x + 1))
    log_expression_df.insert(loc=0, column='sampleid', value=sample_ids)

    # Calculate expression residuals
    def calculate_residuals(gene_id):
        """Calculate gene residuals"""
        gene = gene_id
        exprs_val = log_expression_df[['sampleid', gene]]
        test_df = exprs_val.merge(covariate_df, on='sampleid', how='left')
        test_df = test_df.rename(columns={test_df.columns[1]: 'expression'})
        y, x = dmatrices(
            'expression ~ sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + age + pf1 + pf2',
            test_df,
        )
        model = sm.OLS(y, x)
        residuals = list(model.fit().resid)
        return residuals

    residual_df = pd.DataFrame(list(map(calculate_residuals, gene_ids))).T
    residual_df.columns = gene_ids
    residual_df = residual_df.assign(sampleid=list(sample_ids))
    # residual_df.to_csv(f'gs://{OUTPUT_BUCKET}/kat/chr22_log_residuals.tsv')

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

    gene_info = geneloc_df.iloc[idx]
    snps_within_region = snploc_df[
        snploc_df['pos'].between(gene_info['left'], gene_info['right'])
    ]
    gene_snp_df = snploc_df.merge(pd.DataFrame(snps_within_region))
    gene_snp_df = gene_snp_df.assign(geneid=gene_info.geneid)
    spearman_df = pd.DataFrame(list(gene_snp_df.apply(spearman_correlation, axis=1)))
    spearman_df.columns = ['geneid', 'snpid', 'coef', 'p.value']
    # add in global position and round
    locus = spearman_df.snpid.str.split('_', expand=True)[0]
    chromosome = 'chr' + locus.str.split(':', expand=True)[0]
    position = locus.str.split(':', expand=True)[1]
    spearman_df['locus'], spearman_df['chromosome'], spearman_df['position'] = [
        locus,
        chromosome,
        position,
    ]
    spearman_df['round'] = 1
    # convert to hail table
    # spearman_df.to_csv(f'gs://{OUTPUT_BUCKET}/kat/spearman_df.csv')
    # hl.init_local(default_reference='GRCh38')
    # t = hl.import_table(f'gs://{OUTPUT_BUCKET}/kat/spearman_df.csv', delimiter=',', types={'position': hl.tint32, 'coef': hl.tfloat64, 'p.value': hl.tfloat64}) # noqa: E501; pylint: disable=line-too-long
    # t = t.annotate(global_position=hl.locus(t.chromosome, t.position).global_position()) # noqa: E501; pylint: disable=line-too-long
    # get alleles
    # mt.rows()[c.liftover].alleles
    # turn back into pandas df
    # spearman_df = t.to_pandas()
    return spearman_df


backend = hb.ServiceBackend(billing_project='tob-wgs', bucket='cpg-tob-wgs-test')
b = hb.Batch(name='eQTL', backend=backend, default_python_image=DRIVER_IMAGE)

spearman_dfs_from_scatter = []
# for i in range(get_number_of_scatters()):
for i in range(5):
    j = b.new_python_job(name=f'process_{i}')
    result: hb.resource.PythonResult = j.call(run_computation_in_scatter, i)
    spearman_dfs_from_scatter.append(result)


def function_that_merges_dataframes(*df_list):
    """Merge all Spearman dfs"""
    print(df_list)
    print(type(df_list))
    merged_df: pd.DataFrame = pd.concat(df_list)
    pvalues = merged_df['p.value']
    fdr_values = pd.DataFrame(list(multi.fdrcorrection(pvalues))).iloc[1]
    merged_df = merged_df.assign(FDR=fdr_values)
    merged_df['FDR'] = merged_df.FDR.astype(float)
    return merged_df.to_string()


merge_job = b.new_python_job(name='merge_scatters')
result_second = merge_job.call(
    function_that_merges_dataframes, *spearman_dfs_from_scatter
)
b.write_output(result_second.as_str(), f'gs://{OUTPUT_BUCKET}/kat/test_log.csv')
b.run()
