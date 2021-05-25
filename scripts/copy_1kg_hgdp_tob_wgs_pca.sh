#!/bin/bash

set -ex

# Copy a subset of the data to the `test` bucket.
gsutil -m cp -r gs://cpg-tob-wgs-analysis/1kg_hgdp_tobwgs_pca/v1/hgdp1kg_tobwgs_joined_all_samples.mt gs://cpg-tob-wgs-test/1kg_hgdp_tobwgs_pca/v1/
gsutil -m cp -r gs://cpg-tob-wgs-analysis/1kg_hgdp_nfe/v0/loadings.ht gs://cpg-tob-wgs-analysis/1kg_hgdp_nfe/v0/
