# ancestry

Prioritising and assigning ancestry to populations.

## Running hail

To connect to hail with the required gnomAD dependencies, run the following command:

```hailctl dataproc start \
--max-age 2h \
--region australia-southeast1 \
--packages gnomad \
my-cluster```

Some commands launch jobs that are computationally expensive and require more workers. For these jobs, you can modify the cluster. I specify 100 preemptible workers with the following command:

`hailctl dataproc modify my-cluster --num-preemptible-workers 100 --num-workers 2`

To stop the cluster, run the following:

```hailctl dataproc stop```

## Using and viewing Jupyter notebooks

To launch a Jupyter notebook, run the following command:

```hailctl dataproc connect my-cluster nb```

GitHub renders Jupyter notebooks as static HTML files, so _interactive_ plots [cannot](https://docs.github.com/en/github/managing-files-in-a-repository/working-with-jupyter-notebook-files-on-github) be viewed directly through the repository. In order to view all the interactive features of the notebook, you can enter the location of the `.ipynb` file in the [nbviewer](https://nbviewer.jupyter.org/).

To view the HGDP_1KG_HailExploration.ipynb Jupyter notebook, [click here.](https://nbviewer.jupyter.org/github/populationgenomics/ancestry/blob/main/scripts/HGDP_1KG_HailExploration.ipynb)
