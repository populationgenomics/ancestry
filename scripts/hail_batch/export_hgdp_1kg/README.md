# Export HGDP/1kG subset to PLINK format

This runs a Hail query script in Dataproc using Hail Batch in order to export the a subset of the HGDP/1kG dataset to PLINK format. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level test --output-dir "1kg_hgdp_plink/v0" \
--description "1kg-hgdp-plink" python3 main.py
```
