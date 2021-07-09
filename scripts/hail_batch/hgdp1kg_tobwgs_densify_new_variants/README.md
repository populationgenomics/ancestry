# Densify TOB-WGS data before filtering for new variants

This runs a Hail query script in Dataproc using Hail Batch in order to densify the TOB-WGS matrix before filtering for variants which were selected specifically for the combined TOB-WGS + HGDP/1KG datasets (using the same filtering criteria as gnomAD v2.1, however not limited to exonic regions). To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "1kg_hgdp_densify_new_variants/v0" \
--description "densify tob-wgs new variants" python3 main.py
```
