# Plot PCA of combined TOB-WGS + SNP-chip data

This runs a Hail query script in Dataproc using Hail Batch in order to plot a PCA of the combined TOB-WGS/SNP-chip datasets, generated from the script `snp_chip_generate_pca/snp_chip_generate_pca.py`. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "tob_wgs_snp_chip_pca/v0" \
--description "pca_tob_wgs_snp_chip" python3 main.py
```
