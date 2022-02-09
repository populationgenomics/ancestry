#!/usr/bin/env python3

"""Run Spearman rank correlation on SNPs and expression residuals"""

import os
import hail as hl
import hailtop.batch as hb
import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.stats.multitest as multi
from patsy import dmatrices  # pylint: disable=no-name-in-module
from scipy.stats import spearmanr
import click

DEFAULT_DRIVER_MEMORY = '4G'
DEFAULT_DRIVER_IMAGE = 'australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d2a9c316d6d752edb27623542c8a062db4466842-hail-0.2.73.devc6f6f09cec08'  # noqa: E501; pylint: disable=line-too-long
DRIVER_IMAGE = os.getenv('DRIVER_IMAGE', DEFAULT_DRIVER_IMAGE)


def filter_lowly_expressed_genes(expression_df):
    """Remove genes with low expression in all samples"""

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

    return expression_df


def get_number_of_scatters(expression_df, geneloc_df):
    """get index of total number of genes"""

    expression_df = filter_lowly_expressed_genes(expression_df)
    gene_ids = list(expression_df.columns.values)[1:]
    geneloc_df = geneloc_df[geneloc_df.geneid.isin(gene_ids)]

    return len(geneloc_df.index)


def get_log_expression(expression_df):

    """get logged expression values"""

    expression_df = filter_lowly_expressed_genes(expression_df)

    # Prepare variables
    sample_ids = expression_df.iloc[:, 0]

    # log expression values
    to_log = expression_df.iloc[:, 1:].columns
    log_expression_df = expression_df[to_log].applymap(lambda x: np.log(x + 1))
    log_expression_df.insert(loc=0, column='sampleid', value=sample_ids)

    return log_expression_df


def calculate_residuals(expression_df, covariate_df, output_prefix):
    """Calculate residuals for each gene in scatter"""

    log_expression_df = get_log_expression(expression_df)
    # Prepare variables
    gene_ids = list(log_expression_df.columns.values)[1:]
    sample_ids = log_expression_df.iloc[:, 0]

    # Calculate expression residuals
    def calculate_gene_residual(gene_id):
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

    residual_df = pd.DataFrame(list(map(calculate_gene_residual, gene_ids))).T
    residual_df.columns = gene_ids
    residual_df = residual_df.assign(sampleid=list(sample_ids))
    residual_df.to_csv(os.path.join(output_prefix, f'log_residuals.tsv'))

    return residual_df


# Run Spearman rank in parallel by sending genes in a batches
def run_spearman_correlation_scatter(
    idx,
    expression_df,
    genotype_df,
    geneloc_df,
    snploc_df,
    residuals_df,
):  # pylint: disable=too-many-locals
    """Run genes in scatter"""

    log_expression_df = get_log_expression(expression_df)

    # Prepare variables used to calculate Spearman's correlation
    gene_ids = list(log_expression_df.columns.values)[1:]
    genotype_df = genotype_df[genotype_df.sampleid.isin(log_expression_df.sampleid)]

    # Get 1Mb sliding window around each gene
    geneloc_df = geneloc_df[geneloc_df.geneid.isin(gene_ids)]
    geneloc_df = geneloc_df.assign(left=geneloc_df.start - 1000000)
    geneloc_df = geneloc_df.assign(right=geneloc_df.end + 1000000)
    # geneloc_df.to_csv(output_prefix + f'_gene_SNP_pairs.tsv')

    def spearman_correlation(df):
        """get Spearman rank correlation"""
        gene = df.geneid
        snp = df.snpid
        res_val = residuals_df[['sampleid', gene]]
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
    chromosome = locus.str.split(':', expand=True)[0]
    position = locus.str.split(':', expand=True)[1]
    spearman_df['locus'], spearman_df['chromosome'], spearman_df['position'] = [
        locus,
        chromosome,
        position,
    ]
    spearman_df['round'] = 1
    # convert to hail table. Can't call `hl.from_pandas(spearman_df)` directly
    # because it doesnt' work with the spark local backend
    spearman_df.to_csv(f'spearman_df_{idx}.csv')
    hl.init(default_reference='GRCh37')
    t = hl.import_table(
        f'spearman_df_{idx}.csv',
        delimiter=',',
        types={'position': hl.tint32, 'coef': hl.tfloat64, 'p.value': hl.tfloat64},
    )
    t = t.annotate(global_position=hl.locus(t.chromosome, t.position).global_position())
    # get alleles
    # mt.rows()[c.liftover].alleles
    # turn back into pandas df. Can't call `spearman_df = t.to_pandas()` directly
    # because it doesn't work with the spark local backend
    t.export(f'spearman_df_annotated_{idx}.tsv')
    spearman_df = pd.read_csv(f'spearman_df_annotated_{idx}.tsv', sep='\t')
    return spearman_df


def merge_df_and_convert_to_string(*df_list):
    """Merge all Spearman dfs and convert to string using .to_string() on df"""
    merged_df: pd.DataFrame = pd.concat(df_list)
    pvalues = merged_df['p.value']
    fdr_values = pd.DataFrame(list(multi.fdrcorrection(pvalues))).iloc[1]
    merged_df = merged_df.assign(FDR=fdr_values)
    merged_df['FDR'] = merged_df.FDR.astype(float)
    return merged_df.to_string()


# Create click command line to enter dependency files
@click.command()
@click.option(
    '--expression', required=True, help='A sample x gene TSV of expression values'
)
@click.option('--genotype', required=True, help='A TSV of genotypes for each sample')
@click.option(
    '--geneloc', required=True, help='A TSV of start and end positions for each gene'
)
@click.option(
    '--snploc',
    required=True,
    help='A TSV of snp IDs with chromsome and position values for each',
)
@click.option(
    '--covariates', required=True, help='A TSV of covariates to calculate residuals'
)
@click.option(  # pylint: disable=too-many-locals
    '--output-prefix',
    required=True,
    help='A path prefix of where to output files, eg: gs://MyBucket/output-folder/',
)
def main(
    expression: str,
    genotype,
    geneloc,
    snploc,
    covariates,
    output_prefix: str,
):
    """
    Creates a Hail Batch pipeline for calculating EQTLs
    """
    dataset = os.getenv('DATASET')
    access_level = os.getenv('ACCESS_LEVEL')
    backend = hb.ServiceBackend(
        billing_project=dataset, bucket=f'cpg-{dataset}-{access_level}'
    )
    batch = hb.Batch(name='eQTL', backend=backend, default_python_image=DRIVER_IMAGE)

    # load in files literally to do the get_number of scatters
    expression_df_literal = pd.read_csv(expression, sep='\t')
    geneloc_df_literal = pd.read_csv(geneloc, sep='\t')

    # load files into a python job to avoid memory issues during a submission
    load_job = batch.new_python_job('load-data')
    expression_df = load_job.call(pd.read_csv, expression, sep='\t')
    genotype_df = load_job.call(pd.read_csv, genotype, sep='\t')
    geneloc_df = load_job.call(pd.read_csv, geneloc, sep='\t')
    snploc_df = load_job.call(pd.read_csv, snploc, sep='\t')
    covariate_df = load_job.call(pd.read_csv, covariates, sep='\t')
    calculate_residuals_job = batch.new_python_job('load-data')
    residuals_df = calculate_residuals_job.call(
        calculate_residuals,
        expression_df=expression_df,
        covariate_df=covariate_df,
        output_prefix=output_prefix,
    )

    spearman_dfs_from_scatter = []
    for idx in range(get_number_of_scatters(expression_df_literal, geneloc_df_literal)):
        # for i in range(5):
        j = batch.new_python_job(name=f'process_{idx}')
        result: hb.resource.PythonResult = j.call(
            run_spearman_correlation_scatter,
            idx=idx,
            expression_df=expression_df,
            genotype_df=genotype_df,
            geneloc_df=geneloc_df,
            snploc_df=snploc_df,
            residuals_df=residuals_df,
        )
        spearman_dfs_from_scatter.append(result)

    merge_job = batch.new_python_job(name='merge_scatters')
    result_second = merge_job.call(
        merge_df_and_convert_to_string, *spearman_dfs_from_scatter
    )
    corr_result_output_path = os.path.join(output_prefix + 'correlation_results.csv')
    batch.write_output(result_second.as_str(), corr_result_output_path)
    batch.run(wait=False)


if __name__ == '__main__':
    # pylint: disable=no-value-for-parameter
    main()
