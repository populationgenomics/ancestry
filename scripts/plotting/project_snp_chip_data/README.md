# Project SNP-chip data onto HGDP/1kG samples

This runs a Hail query script in Dataproc using Hail Batch in order project the SNP-Chip samples onto the HGDP/1kG dataset. To run, use conda to install the analysis-runner, then execute the following command:

```sh
POP=nfe
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "gs://cpg-tob-wgs-main-web/1kg_hgdp_${POP}_snp_chip/v0" \
--description "${POP} snp-chip" python3 main.py ${POP}
```
