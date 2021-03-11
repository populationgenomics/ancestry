# ancestry

Prioritising and assigning ancestry to populations.

## Running hail

To connect  to hail with the required gnomAD dependencies, run the following command:

```hailctl dataproc start --max-age 2h --region australia-southeast1 my-cluster --packages gnomad```

Some commands are computationally expensive and require more workers. For these jobs, you can modify the cluster. I specify 100 preemptible workers with the following command:

```hailctl dataproc modify my-cluster --num-preemptible-workers 100 --num-workers 2```

To stop the cluster, run the following:

```hailctl dataproc stop```

## Using and viewing Jupyter notebooks

To launch a Jupyter notebook, run the following command:

```hailctl dataproc connect my-cluster nb```

Not all images in Jupyter notebooks are able to be viewed on Github. In order to view all of the Jupyter notebook images, use nbviewer instead of clicking on the Jupyter notebook file within this Github repository.

To view the script HGDP_1KG_HailExploration.ipynb, (click here.)[https://nbviewer.jupyter.org/github/populationgenomics/ancestry/blob/main/scripts/HGDP_1KG_HailExploration.ipynb]
