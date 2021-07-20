# Project WGS samples onto SNP-chip PCA

This runs a Hail query script in Dataproc using Hail Batch in order to generate a PCA of the TOB SNP-chip data, then project the WGS samples onto it. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "tob_wgs_snp_chip_pca/v2" \
--description "project wgs smaples" python3 main.py
```
