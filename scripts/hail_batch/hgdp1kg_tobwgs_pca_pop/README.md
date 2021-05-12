# Perform pca and variant analysis on pop-specific samples

This runs a Hail query script in Dataproc using Hail Batch in order to perform pca and variant frequency analysis on a pop-speific samples (e.g. non-Finnish European (nfe)) from the 1KG + HGDP dataset. To run, use conda to install the analysis-runner, then execute the following command:

```sh
POP=nfe
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "gs://cpg-tob-wgs-analysis/1kg_hgdp_${POP}/v0" \
--description "pca ${POP}" python3 main.py ${POP}
```

Depends on `hgdp1kg_tobwgs_pca/hgdp_1kg_tob_wgs_pca.py`.
