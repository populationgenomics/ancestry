# Generate PCA for NFE samples

This runs a Hail query script in Dataproc using Hail Batch in order to generate a PCA for NFE samples from the final TOB batch (PBMC + bone marrow samples), with all outlier samples removed. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "tob_wgs_pca/nfe_no_outliers/v0" \
--description "tob nfe pca" python3 main.py
```
