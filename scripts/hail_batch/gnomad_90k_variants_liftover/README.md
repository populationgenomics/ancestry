# Perform liftover on the gnomad v2 loadings

This runs a Hail query script in Dataproc using Hail Batch in order to perform liftover on the gnomAD v2 loadings (GRCh37 -> GRCh38). To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset ancestry \
--access-level test --output-dir "gs://cpg-reference-temporary/gnomad_loadings_liftover/v0" \
--description "gnomad liftover" main.py
```
