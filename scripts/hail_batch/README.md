# Generate PCA loadings, eigenvalues, and scores
This runs a Hail query script in Dataproc using Hail Batch in order to output the PCA eigenvalues, scores, and loadings from the 1KG + HGDP dataset. To run, use conda to install the analysis-runner, then execute the following command:

```
analysis-runner --dataset ancestry \
--access-level test --output-dir “gs://cpg-ancestry-main/Kat_PCA_loadings" \
--description “PCA loadings” main.py
```


