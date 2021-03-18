#!/usr/bin/env python

# # Hail PCA Exploration on HGDP + 1KG data

from bokeh.models import CategoricalColorMapper
from bokeh.palettes import turbo
from pprint import pprint
from gnomad.utils.annotations import annotate_adj
from hail.plot import show

import hail as hl
hl.init()
hl.plot.output_notebook()

# read in table
mt = hl.read_matrix_table("gs://gcp-public-data--gnomad/release/3.1/mt/genomes/gnomad.genomes.v3.1.hgdp_1kg_subset_dense.mt")
# data epxloration
mt = mt.head(100, n_cols = 100)
mt.count_rows()
mt.count_cols()
mt.cols().show()
mt.rows().show()

# perform sampole QC
mt_qc = hl.sample_qc(mt)
p = hl.plot.histogram(mt_qc.sample_qc.call_rate, range=(.88,1), legend='Call Rate')
p_2 = hl.plot.histogram(mt_qc.sample_qc.gq_stats.mean, legend='Mean Sample GQ')

# PCA 
columns = mt.cols()
pca_scores = columns.population_inference.pca_scores
labels = columns.population_inference.pop
pops = list(set(labels.collect()))
mapper = CategoricalColorMapper(palette=turbo(8), factors=pops)

# plot each dimension
p1 = hl.plot.scatter(pca_scores[0], pca_scores[1], label=labels,
                    title='PCA', xlabel='PC1', ylabel='PC2', collect_all=True, colors = mapper)

p2 = hl.plot.scatter(pca_scores[1], pca_scores[2], label=labels,
                    title='PCA', xlabel='PC2', ylabel='PC3', collect_all=True, colors = mapper)

p3 = hl.plot.scatter(pca_scores[2], pca_scores[3], label=labels,
                    title='PCA', xlabel='PC3', ylabel='PC4', collect_all=True, colors = mapper)


p4 = hl.plot.scatter(pca_scores[3], pca_scores[4], label=labels,
                    title='PCA', xlabel='PC4', ylabel='PC5', collect_all=True, colors = mapper)

p5 = hl.plot.scatter(pca_scores[4], pca_scores[5], label=labels,
                    title='PCA', xlabel='PC5', ylabel='PC6', collect_all=True, colors = mapper)
