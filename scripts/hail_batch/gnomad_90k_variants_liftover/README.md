# Perform liftover on the gnomad v2 loadings

This runs a Hail query script in Dataproc using Hail Batch in order to perform liftover on the gnomAD v2 loadings (GRCh37 -> GRCh38). To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset ancestry \
--access-level test --output-dir "gs://cpg-ancestry-temporary/gnomad_loadings_liftover_2.1/v0" \
--description "gnomad liftover 2.1" main.py
```
