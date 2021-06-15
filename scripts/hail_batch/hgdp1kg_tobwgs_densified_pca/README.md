# Densify TOB-WGS data

This runs a Hail query script in Dataproc using Hail Batch in order to perform PCA on the densified TOB-WGS matrix table. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "gs://cpg-tob-wgs-main/1kg_hgdp_densified_pca/v0" \
--description "PCA on densified tob-wgs" python3 main.py
```
