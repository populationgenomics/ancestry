# Generate PCA for NFE samples

This runs a Hail query script in Dataproc using Hail Batch in order to generate a PCA for NFE samples from the final TOB batch, with all outlier samples removed. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level test --output-dir "kat/pca/nfe/v0" \
--description "tob nfe pca" python3 main.py
```
