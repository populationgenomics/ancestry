# Densify benchmark

This runs a Hail query script in Dataproc using Hail Batch in order to test the cost and speed of densifying 20 samples. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "densify_benchmark/v0" \
--description "densify benchmark" python3 main.py
```
