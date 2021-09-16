# Repartition the lifted over gnomAD v2 loadings

This runs a Hail query script in Dataproc using Hail Batch in order to repartition the lifted over gnomAD v2 loadings, generated in `gnomad_90k_variants_liftover/gnomad_loadings_90k_liftover.py`. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "gnomad/v0" \
--description "repartition-loadings" python3 main.py
``
