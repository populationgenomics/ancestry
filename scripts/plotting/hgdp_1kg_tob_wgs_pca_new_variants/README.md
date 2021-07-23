# Plot PCA and loadings of newly-selected variants

This runs a Hail query script in Dataproc using Hail Batch in order to make a PCA plot of the HGDP/1kG + TOB-WGS dataset, filtered by the variants selected from the script `variant_selection/hgdp_1kg_tob_wgs_variant_selection.py`. This file is dependent on the PCA generation step, located in `hgdp1kg_tobwgs_densified_pca_new_variants/hgdp_1kg_tob_wgs_densified_pca_new_variants.py`. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "1kg_hgdp_densified_pca_new_variants/v0" \
--description "plot pca new variants" python3 main.py
```
