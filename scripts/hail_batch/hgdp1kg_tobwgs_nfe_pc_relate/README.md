# Run `pc_relate` on nfe samples

This runs a Hail query script in Dataproc using Hail Batch in order to estimate kinship on nfe samples using the `pc_relate` function in Hail. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "tob_wgs_hgdp_1kg_nfe_pc_relate/v0" \
--description "nfe pc_relate" python3 main.py
```
