# Perform pca on the three 1kG reprocessed GATK4 samples

This runs a Hail query script in Dataproc using Hail Batch in order to plot a PCA of the three 1kG samples reprocessed using the GATK4 warp pipeline. To run, use conda to install the analysis-runner, then execute the following command:

```sh
POP=nfe
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "gs://cpg-tob-wgs-web/1kg_hgdp_${POP}_reprocessed_warp/v0" \
--description "pca ${POP} reprocessed warp" python3 main.py ${POP}
```
