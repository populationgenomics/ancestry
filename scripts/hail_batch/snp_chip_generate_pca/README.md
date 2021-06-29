# Generate PCA of variants in SNP-chip data

This runs a Hail query script in Dataproc using Hail Batch in order to generate a PCA of the TOB-WGS data, but filtered to variants present in the TOB SNP-chip data. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "gs://cpg-tob-wgs-main/tob_wgs_snp_chip_variant_pca/v0" \
--description "snp chip variants pca" python3 main.py
```
