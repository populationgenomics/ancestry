# Run hail batch job to generate list of eQTLs

This runs a Hail batch script in order to generate a list of eQTLs from scRNA-seq expression. This code was taken from Seyhan Yazar from Joseph Powell's group at the Garvan-Weizmann Centre for Cellular Genomics, then converted into Python/hail batch. To run, use conda to install the analysis-runner. For the first round of eQTLs, execute the following command:

```sh
analysis-runner --dataset tob-wgs \
    --access-level test --output-dir "kat/v0" \
    --description "eqtl batch job" \
    python3 generate_eqtl_spearman.py \
        --output-prefix 'gs://cpg-tob-wgs-test/kat/plasma/chr22/'
        --expression 'gs://cpg-tob-wgs-test/kat/input/Plasma_expression.tsv' \
        --genotype 'gs://cpg-tob-wgs-test/kat/input/genotype_chr22.tsv' \
        --geneloc 'gs://cpg-tob-wgs-test/kat/input/geneloc_chr22.tsv' \
        --snploc 'gs://cpg-tob-wgs-test/kat/input/snpsloc_chr22.tsv' \
        --covariates 'gs://cpg-tob-wgs-test/kat/input/Plasma_peer_factors.tsv' \
```

For the conditional analysis (rounds 2-5), execute the following command:

```sh
analysis-runner --dataset tob-wgs \
    --access-level test --output-dir "kat/v0" \
    --description "eqtl batch job" \
    python3 round2.conditional_analysis_test.py \
        --output_prefix 'gs://cpg-tob-wgs-test/kat/plasma/chr22/' \
        --residuals 'gs://cpg-tob-wgs-test/kat/plasma/chr22/log_residuals.tsv' \
        --significant_snps 'gs://cpg-tob-wgs-test/kat/plasma/chr22/correlation_results.csv' \
        --genotype 'gs://cpg-tob-wgs-test/kat/input/genotype_chr22.tsv' \
        --test_subset_genes 5 # test with 5 genes only
```
