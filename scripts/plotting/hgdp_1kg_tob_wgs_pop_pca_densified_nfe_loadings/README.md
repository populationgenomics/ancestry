# Plot PCA loadings of densified PCA on NFE samples

This runs a Hail query script in Dataproc using Hail Batch in order to plot the output of the PCA loadings on nfe samples only, generated from the script `hgdp_1kg_tob_wgs_densified_pca.py`. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "gs://cpg-tob-wgs-main-web/tob_wgs_hgdp_1kg_nfe_pca_densified/v3" \
--description "densified pca loadings nfe" python3 main.py
```
