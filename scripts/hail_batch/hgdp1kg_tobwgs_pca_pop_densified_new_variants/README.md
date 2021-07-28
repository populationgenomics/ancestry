# Perform pca on pop-specific samples using newly-selected variants

This runs a Hail query script in Dataproc using Hail Batch in order to perform PCA for population subsets on the densified TOB-WGS matrix table, filtered for variants selected specifically for the combined TOB-WGS + HGDP/1KG datasets (using the same filtering criteria as gnomAD v2.1, however not limited to exonic regions). To run, use conda to install the analysis-runner, then execute the following command:

```sh
POP=nfe
analysis-runner --dataset tob-wgs \
--access-level standard --output-dir "1kg_hgdp_densified_${POP}_new_variants/v0" \
--description "pca ${POP} densified" python3 main.py ${POP}
```

Depends on `hgdp1kg_tobwgs_densified_pca_new_variants/hgdp_1kg_tob_wgs_densified_pca_new_variants.py`.
