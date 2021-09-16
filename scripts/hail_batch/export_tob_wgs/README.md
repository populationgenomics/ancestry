# Export TOB-WGS dataset to PLINK format

This runs a Hail query script in Dataproc using Hail Batch in order to export the TOB-WGS joint callset to PLINK format. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "tob_wgs_plink/v0" \
--description "TOB plink" python3 main.py
```
