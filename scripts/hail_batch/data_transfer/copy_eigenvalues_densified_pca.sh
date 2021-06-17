#!/usr/bin/env bash

set -ex

gsutil cp "gs://cpg-tob-wgs-main/1kg_hgdp_densified_pca/v2/eigenvalues.ht" $OUTPUT
