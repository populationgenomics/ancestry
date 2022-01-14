# Calculate LD for each locus within the TOB dataset 

This runs a Hail query script in Dataproc using Hail Batch in order to calculate the LD score for each locus within the TOB dataset. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level test --output-dir "kat/v0" \
--description "ld-calculate" python3 main.py
```
