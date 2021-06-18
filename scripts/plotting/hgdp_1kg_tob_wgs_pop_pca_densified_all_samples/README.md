# Plot PCA of densified PCA on all HGDP/1kG + TOB-WGS samples

This runs a Hail query script in Dataproc using Hail Batch in order to plot the output of the HGDP/1kG + TOB-WGS PCA results, generated from the script `hgdp_1kg_tob_wgs_densified_pca.py`. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "gs://cpg-tob-wgs-main-web/tob_wgs_hgdp_1kg_pca_densified/v0" \
--description "densified pca all samples" python3 main.py
```
