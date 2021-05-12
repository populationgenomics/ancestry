# Plotting

First, setup the environment:

```sh
conda env create -n ancestry --file environment.yml
conda activate ancestry
# To allow hail reading from GCS buckets:
curl -sSL https://broad.io/install-gcs-connector | python3
```

Run the script locally:

```sh
python hgdp_1kg_tob_wgs_plot_pca.py --number-of-pcs 5
```
