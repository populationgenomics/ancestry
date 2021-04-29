# Perform a pca on the combined HGDP + 1KG and TOB-WGS data

This runs a Hail query script in Dataproc using Hail Batch in order to perform a pca on the 1KG + HGDP + TOB-WGS dataset. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "gs://cpg-tob-wgs-analysis/1kg_hgdp_tobwgs_pca/v0" \
--description "hgdp1kg tobwgs pca" python3 main.py
```
