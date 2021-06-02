# Test the densify function on the TOB-WGS data

This runs a Hail query script in Dataproc using Hail Batch in order to test the densify function on the TOB-WGS dataset. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level test --output-dir "gs://cpg-tob-wgs-test/densify/v0" \
--description "tobwgs densify test" python3 main.py
```
