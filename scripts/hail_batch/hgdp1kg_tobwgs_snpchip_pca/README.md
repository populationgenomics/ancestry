# Densify TOB-WGS data

This runs a Hail query script in Dataproc using Hail Batch in order to perform PCA on the combined TOB-WGS + HGDP/1kG data. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level test --output-dir "gs://cpg-tob-wgs-test/tob_wgs_snpchip_pca/v0" \
--description "snpchip pca" python3 main.py
```
