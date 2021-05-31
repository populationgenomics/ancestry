# Variant selection exploration for HGDP/1kG + TOB-WGS data

This runs a Hail query script in Dataproc using Hail Batch in order to select a set of variants from the HGDP/1kG + TOB-WGS data, based off of gnomAD v3 parameters. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "gs://cpg-tob-wgs-analysis/tob_wgs_hgdp_1kg_variant_selection_exploration/v0" \
--description "variant selection exploration" python3 main.py
```
