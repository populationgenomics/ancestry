# histogram of newly-selected variants

This runs a Hail query script in Dataproc using Hail Batch in order to make a histogram of the variants selected from the script `variant_selection/hgdp_1kg_tob_wgs_variant_selection.py`. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "tob_wgs_hgdp_1kg_variant_selection/v0" \
--description "variant selection histogram" python3 main.py
```
