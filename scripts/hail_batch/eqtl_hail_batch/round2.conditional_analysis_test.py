"""Perform conditional analysis on SNPs and expression residuals"""

import os

import hail as hl
import hailtop.batch as hb
import pandas as pd
import statsmodels.api as sm
import statsmodels.stats.multitest as multi
from patsy import dmatrices  # pylint: disable=no-name-in-module
from scipy.stats import spearmanr

import click


DEFAULT_DRIVER_MEMORY = '4G'
DEFAULT_DRIVER_IMAGE = 'australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d2a9c316d6d752edb27623542c8a062db4466842-hail-0.2.73.devc6f6f09cec08'  # noqa: E501; pylint: disable=line-too-long
DRIVER_IMAGE = os.getenv('DRIVER_IMAGE', DEFAULT_DRIVER_IMAGE)


def get_number_of_scatters(residual_df, significant_snps_df):
    """get index of total number of genes"""

    # Identify the top eSNP for each eGene and assign remaining to df
    esnp1 = (
        significant_snps_df.sort_values(['geneid', 'p.value'], ascending=True)
        .groupby('geneid')
        .first()
        .reset_index()
    )
    gene_ids = esnp1['geneid'][esnp1['geneid'].isin(residual_df.columns)]

    return len(gene_ids)


def get_genotype_df(genotype_df, significant_snps_df, sample_ids):
    """load genotype df and filter"""
    # Subset genotype file for the significant SNPs
    genotype_df = genotype_df.loc[
        genotype_df['sampleid'].isin(  # pylint: disable=unsubscriptable-object
            sample_ids.sampleid
        ),
        :,
    ]
    genotype_sampleid = genotype_df['sampleid']
    genotype_df = genotype_df.loc[
        :, genotype_df.columns.isin(significant_snps_df.snpid)
    ]
    genotype_df = genotype_df.assign(sampleid=genotype_sampleid)

    return genotype_df


def calculate_residual_df(genotype_df, residual_df, significant_snps_df):
    """calculate residuals for gene list"""

    # Identify the top eSNP for each eGene and assign remaining to df
    esnp1 = (
        significant_snps_df.sort_values(['geneid', 'FDR'], ascending=True)
        .groupby('geneid')
        .first()
        .reset_index()
    )

    # Subset residuals for the genes to be tested
    sample_ids = residual_df.loc[:, ['sampleid']]
    gene_ids = esnp1['geneid'][esnp1['geneid'].isin(residual_df.columns)]
    residual_df = residual_df.loc[:, residual_df.columns.isin(gene_ids)]
    residual_df['sampleid'] = sample_ids
    genotype_df = get_genotype_df(genotype_df, significant_snps_df, sample_ids)

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

    # mfranklin: I've commented this out, because iteration isn't defined
    #             and it's written outside this loop. Maybe this could be removed?
    # adjusted_residual_mat.to_csv(
    #     f'gs://{OUTPUT_BUCKET}/kat/eSNP1_adjusted_residuals_{iteration}.tsv'
    # )

    return adjusted_residual_mat


# Run Spearman rank in parallel by sending genes in a batches
def run_computation_in_scatter(
    iteration,  # pylint: disable=redefined-outer-name
    idx,
    genotype_df,
    residual_df,
    significant_snps_df,
):
    """Run genes in scatter"""
    # Input filenames

    print(f'iteration = {iteration}')
    print(f'idx = {idx}')

    # make sure 'geneid' is the first column
    # otherwise, error thrown when using reset_index
    cols = list(significant_snps_df)
    cols.insert(0, cols.pop(cols.index('geneid')))
    significant_snps_df = significant_snps_df.loc[:, cols]
    esnps_to_test = (
        significant_snps_df.sort_values(['geneid', 'FDR'], ascending=True)
        .groupby('geneid')
        .apply(lambda group: group.iloc[1:, 1:])
        .reset_index()
    )

    sample_ids = residual_df.loc[:, ['sampleid']]
    genotype_df = get_genotype_df(genotype_df, significant_snps_df, sample_ids)

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
        significant_snps_df.sort_values(['geneid', 'FDR'], ascending=True)
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
    locus = adjusted_spearman_df.snpid.str.split('_', expand=True)[0]
    chromosome = locus.str.split(':', expand=True)[0]
    position = locus.str.split(':', expand=True)[1]

    adjusted_spearman_df['locus'] = locus
    adjusted_spearman_df['chromosome'] = chromosome
    adjusted_spearman_df['position'] = position
    adjusted_spearman_df['round'] = iteration

    # convert to hail table. Can't call `hl.from_pandas(spearman_df)` directly
    # because it doesnt' work with the spark local backend
    adjusted_spearman_df.to_csv('adjusted_spearman_df.csv')
    hl.init(default_reference='GRCh37')
    t = hl.import_table(
        'adjusted_spearman_df.csv',
        delimiter=',',
        types={'position': hl.tint32, 'coef': hl.tfloat64, 'p.value': hl.tfloat64},
    )
    t = t.annotate(
        global_position=hl.locus(t.chromosome, t.position).global_position()
    )  # noqa: E501; pylint: disable=line-too-long
    # get alleles and rsid
    # mt = hl.read_matrix_table(HGDP1KG_TOBWGS)
    t = t.annotate(locus=hl.locus(t.chromosome, t.position))
    # mt.rows()[t.locus].alleles
    # mt.rows()[t.locus].rsid
    t = t.annotate(
        id=hl.str(':').join(
            [
                hl.str(t.chromosome),
                hl.str(t.position),
                # result.A1,
                # result.A2,
                t.geneid,
                # result.db_key, # cell_type_id (eg nk, mononc)
                hl.str(t.round),
            ]
        )
    )
    # turn back into pandas df. Can't call `spearman_df = t.to_pandas()` directly
    # because it doesn't work with the spark local backend
    t.export('adjusted_spearman_df_annotated.tsv')
    adjusted_spearman_df = pd.read_csv('adjusted_spearman_df_annotated.tsv', sep='\t')

    # set variables for next iteration of loop
    significant_snps_df = adjusted_spearman_df

    return significant_snps_df


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


def convert_dataframe_to_text(dataframe):
    """
    convert to string for writing
    """
    return dataframe.to_string()


@click.command()
@click.option(
    '--significant_snps', required=True, help='A space separated list of SNPs'
)
@click.option('--residuals', required=True, help='A TSV of gene residuals')
@click.option('--genotype', required=True, help='A TSV of genotypes for each sample')
@click.option(
    '--output_prefix',
    required=True,
    help='A path prefix of where to output files, eg: gs://MyBucket/output-folder/',
)
@click.option(
    '--iterations', type=int, default=5, help='Number of iterations to perform'
)
@click.option(
    '--test_subset_genes',
    type=int,
    help='Test with {test_subset_genes} genes, often = 5.',
)
def main(
    significant_snps: str,
    residuals,
    genotype,
    output_prefix: str,
    iterations=5,
    test_subset_genes=None,
):
    """
    Creates a Hail Batch pipeline for calculating EQTL using {iterations} iterations,
    scattered across the number of genes
    """
    dataset = os.getenv('DATASET')
    access_level = os.getenv('ACCESS_LEVEL')
    backend = hb.ServiceBackend(
        billing_project=dataset, bucket=f'cpg-{dataset}-{access_level}'
    )
    batch = hb.Batch(name='eQTL', backend=backend, default_python_image=DRIVER_IMAGE)

    if test_subset_genes:
        n_genes = test_subset_genes
    else:
        # load these literally to do the get_number of scatters
        print(f'Loading residuals: {residuals}')
        residual_df_literal = pd.read_csv(residuals, sep='\t')
        print(f'Loading significant_snps: {significant_snps}')
        significant_snps_df_literal = pd.read_csv(
            significant_snps, sep=' ', skipinitialspace=True
        )

        print('Loaded data to prepare workflow')
        # test with 5 genes
        n_genes = test_subset_genes or get_number_of_scatters(
            residual_df_literal, significant_snps_df_literal
        )

    # load these in a python job to avoid memory issues during a submission
    load_job = batch.new_python_job('load-data')
    genotype_df = load_job.call(pd.read_csv, genotype, sep='\t')
    residual_df = load_job.call(pd.read_csv, residuals, sep='\t')
    significant_snps_df = load_job.call(
        pd.read_csv, significant_snps, sep=' ', skipinitialspace=True
    )

    previous_sig_snps_result = significant_snps_df  # pylint: disable=invalid-name
    previous_residual_result = residual_df  # pylint: disable=invalid-name
    for iteration in range(iterations):

        calc_resid_df_job = batch.new_python_job(f'calculate-resid-df-iter-{iteration}')
        previous_residual_result = calc_resid_df_job.call(
            calculate_residual_df,
            genotype_df,
            previous_residual_result,
            previous_sig_snps_result,
        )

        # convert residual df to string for output
        residual_as_str = calc_resid_df_job.call(
            convert_dataframe_to_text, previous_residual_result
        )
        # output residual df for each iteration
        batch.write_output(
            residual_as_str.as_str(),
            os.path.join(output_prefix, f'round{iteration+1}_residual_results.csv'),
        )

        sig_snps_dfs = []
        for gene_idx in range(n_genes):
            j = batch.new_python_job(name=f'process_iter_{iteration}_job_{gene_idx}')
            gene_result: hb.resource.PythonResult = j.call(
                run_computation_in_scatter,
                iteration,
                gene_idx,
                genotype_df,
                previous_residual_result,
                previous_sig_snps_result,
            )
            sig_snps_dfs.append(gene_result)

        merge_job = batch.new_python_job(name='merge_scatters')
        previous_sig_snps_result = merge_job.call(
            merge_significant_snps_dfs, *sig_snps_dfs
        )

        # convert sig snps to string for output
        sig_snps_as_string = merge_job.call(
            convert_dataframe_to_text, previous_sig_snps_result
        )
        # output sig snps for each iteration
        sig_snps_output_path = os.path.join(
            output_prefix, f'round{iteration+1}_significant_correlation_results.csv'
        )
        batch.write_output(sig_snps_as_string.as_str(), sig_snps_output_path)

    batch.run(wait=False)


if __name__ == '__main__':
    # pylint: disable=no-value-for-parameter
    main()
