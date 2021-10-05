# Get number of variants with a MAF > 0.05 in TOB dataset

This runs a Hail query script in Dataproc using Hail Batch in order to get the number of variants with a MAF > 0.05. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "tob_maf/v0" \
--description "TOB maf" python3 main.py
```
