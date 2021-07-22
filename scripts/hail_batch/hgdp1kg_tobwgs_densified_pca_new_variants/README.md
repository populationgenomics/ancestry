# PCA on TOB-WGS using newly-selected variants

This runs a Hail query script in Dataproc using Hail Batch in order to perform PCA on the densified TOB-WGS matrix table, filtered for variants selected specifically for the combined TOB-WGS + HGDP/1KG datasets (using the same filtering criteria as gnomAD v2.1, however not limited to exonic regions). To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "1kg_hgdp_densified_pca_new_variants/v0" \
--description "PCA on new tob variants" python3 main.py
```
