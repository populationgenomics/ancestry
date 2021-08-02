# Plot PCA of newly-selected variants on NFE samples without outliers

This runs a Hail query script in Dataproc using Hail Batch in order to plot the output of the PCA results on nfe samples only, generated from the script `hgdp1kg_tobwgs_pca_nfe_new_variants_no_outliers/hgdp_1kg_tob_wgs_nfe_pca_no_outliers.py`. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "tob_wgs_hgdp_1kg_nfe_pca_new_variants/v1" \
--description "pca nfe no outliers" python3 main.py
```
