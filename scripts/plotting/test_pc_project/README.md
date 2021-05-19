# Perform pca on the three 1kG reprocessed GATK4 samples

This runs a Hail query script in Dataproc using Hail Batch in order to test whether the `pc_project` function produces the same results on three samples which were not reprocessed using the GATK4 pipeline as the three samples that were. To run, use conda to install the analysis-runner, then execute the following command:

```sh
POP=nfe
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "gs://cpg-tob-wgs-web/1kg_hgdp_${POP}/v0" \
--description "pc project ${POP}" python3 main.py ${POP}
```
