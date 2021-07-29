# Plot PCA of newly-selected variants on NFE samples

This runs a Hail query script in Dataproc using Hail Batch in order to plot the output of the PCA results on nfe samples only, generated from the script `hgdp1kg_tobwgs_pca_pop_densified_new_variants/hgdp_1kg_tob_wgs_pop_pca_densified.py`. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "tob_wgs_hgdp_1kg_nfe_pca_new_variants/v0" \
--description "pca on nfe" python3 main.py
```
