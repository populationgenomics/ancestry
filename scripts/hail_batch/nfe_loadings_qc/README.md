# Loadings dataset exploration

This runs a Hail query script in Dataproc using Hail Batch in order to explore the nfe loadings. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "gs://cpg-tob-wgs-main-tmp/nfe_loadings_exploration/v0" \
--description "nfe loadings exploration" python3 main.py
```
