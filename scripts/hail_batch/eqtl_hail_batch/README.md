# Run hail batch job to generate list of eQTLs

This runs a Hail batch script in order to generate a list of eQTLs from scRNA-seq expression. This code was taken from Seyhan Yazar from Joseph Powell's group at the Garvan-Weizmann Centre for Cellular Genomics, then converted into Python/hail batch. To run, use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner --dataset tob-wgs \
--access-level test --output-dir "kat/v0" \
--description "eqtl batch job" python3 generate_eqtl_spearman.py
```
