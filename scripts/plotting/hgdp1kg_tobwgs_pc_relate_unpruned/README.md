# Plot related samples from `pc_relate`

This runs a Hail query script in Dataproc using Hail Batch in order to plot kinship estimates after running `pc_relate` on the HGDP/1kG + TOB-WGS samples. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "tob_wgs_hgdp_1kg_pc_relate/v0" \
--description "relaated-samples-plot" python3 main.py
```
