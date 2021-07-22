# Project WGS samples onto SNP-chip PCA

This runs a Hail query script in Dataproc using Hail Batch in order to increase the number of partitions in the SNP-chip matrix table. Since the matrix table currently only has one partition, this results in a high amount of work done using one worker only. This is particularly difficult when intersecting with another matrix table that has many partitions. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "tob_wgs_snp_chip_pca/increase_partitions/v1" \
--description "increase partitions snp-chip" python3 main.py
```
