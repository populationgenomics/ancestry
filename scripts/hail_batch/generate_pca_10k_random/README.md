# Generate PCA loadings, eigenvalues, and scores on 10k randomly-sampled rows

This runs a Hail query script in Dataproc using Hail Batch in order to output the PCA eigenvalues, scores, and loadings from the 1KG + HGDP dataset. The PCA was generated using 10k randomly-selected rows from the original 1KG + HGDP dataset. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset ancestry \
--access-level test --output-dir "gs://cpg-ancestry-temporary/1kg_hgdp_pca/v0" \
--description "PCA loadings 10k" main.py
```
