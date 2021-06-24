# Perform pca on the three 1kG reprocessed GATK4 samples

This runs a Hail query script in Dataproc using Hail Batch in order to plot a PCA of samples reprocessed using the KCCG pipeline. To run, use conda to install the analysis-runner, then execute the following command:

```sh
POP=nfe
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "gs://cpg-tob-wgs-main-web/1kg_hgdp_${POP}_reprocessed_kccg/v0" \
--description "pca ${POP} reprocessed kccg" python3 main.py ${POP}
```
