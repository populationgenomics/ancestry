# Perform pca and variant analysis on NFE samples

This runs a Hail query script in Dataproc using Hail Batch in order to perform pca and variant frequency analysis on non-Finnish European (nfe) samples from the 1KG + HGDP dataset. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset ancestry \
--access-level test --output-dir "gs://cpg-ancestry-temporary/1kg_hgdp_nfe/v0" \
--description "pca nfe" python3 main.py
```
