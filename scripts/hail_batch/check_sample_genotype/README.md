# Perform sample QC on 3 selected TOB-WGS samples

This runs a Hail query script in Dataproc using Hail Batch in order to perform sample QC on three samples, TOB1524, TOB1532, and TOB1533. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level test --output-dir "gs://cpg-tob-wgs-temporary/sample_qc/v0" \
--description "get sample qc metrics" python3 main.py
```
