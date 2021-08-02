# Perform pca on nfe samples using newly-selected variants and removing outlier samples

This runs a Hail query script in Dataproc using Hail Batch in order to perform PCA for nfe samples on the densified TOB-WGS matrix table, filtered for variants selected specifically for the combined TOB-WGS + HGDP/1KG datasets (using the same filtering criteria as gnomAD v2.1, however not limited to exonic regions). Two outlier samples, which were generated in `hgdp_1kg_tob_wgs_pop_pca_new_variants_nfe/hgdp_1kg_tob_wgs_plot_pca_nfe.py`, have been removed in order to see more clear clustering between samples. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "1kg_hgdp_densified_nfe_new_variants/v0" \
--description "pca nfe no outliers" python3 main.py
```

Depends on `hgdp1kg_tobwgs_densified_pca_new_variants/hgdp_1kg_tob_wgs_densified_pca_new_variants.py`.
