# Generate PCA of SNP-chip data only

This runs a Hail query script in Dataproc using Hail Batch in order to generate a PCA of the TOB SNP-chip data. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "tob_snp_chip_pca/v0" \
--description "snp chip only pca" python3 main.py
```
