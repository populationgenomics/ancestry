# Plot PCA of SNP-chip data

This runs a Hail query script in Dataproc using Hail Batch in order to plot a PCA of the TOB SNP-chip dataset, generated from the script `snp_chip_pca_sanity_check/snp_chip_only_generate_pca.py`. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "tob_snp_chip_pca/v0" \
--description "pca_tob_snp_chip" python3 main.py
```
