# Densify TOB-WGS data

This runs a Hail query script in Dataproc using Hail Batch in order to densify the TOB-WGS matrix, then perform PCA. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "gs://cpg-tob-wgs-main/1kg_hgdp_densify/v2" \
--description "densify tob-wgs PCA" python3 main.py
```
