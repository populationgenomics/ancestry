# Perform ld-pruning

This runs a Hail query script in Dataproc using Hail Batch in order to perform ld pruning on the 1KG + HGDP dataset. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset ancestry \
--access-level test --output-dir "gs://cpg-ancestry-temporary/1kg_hgdp_ld_pruning/v0" \
--description "ld pruning" python3 main.py
```
