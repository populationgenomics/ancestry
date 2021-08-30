# Estimate relatedness between NFE samples using King

This runs a Hail query script in Dataproc using Hail Batch in order to estimate relatedness between NFE samples from the HGDP/1KG datasets using the KING relatedness estimator. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "king/v0" \
--description "king nfe samples" python3 main.py
```
