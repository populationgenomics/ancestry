# Densify TOB-WGS data

This runs a Hail query script in Dataproc using Hail Batch in order to densify the TOB-WGS matrix. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level test --output-dir "gs://cpg-tob-wgs-test/1kg_hgdp_densify/v3" \
--description "densify pca test" python3 main.py
```
