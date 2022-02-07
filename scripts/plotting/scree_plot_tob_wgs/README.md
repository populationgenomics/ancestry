# Scree plot of TOB-WGS + NFE samples

This runs a Hail query script in Dataproc using Hail Batch in order to make a scree plot of the TOB-WGS + NFE samples (from HGDP + 1kg dataset), in order to identify the number of PCs to use in the association script LM. This script is dependent on output from `hail_batch/final_tob_batch_pca/generate_pca_no_outliers.py`. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "tob_wgs_hgdp_1kg_nfe/scree/v0" \
--description "scree plot tob" python3 main.py
```
