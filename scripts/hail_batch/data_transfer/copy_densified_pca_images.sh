#!/usr/bin/env bash

set -ex

gsutil cp -r "gs://cpg-tob-wgs-main/tob_wgs_hgdp_1kg_nfe_pca_densified/v0/*.png" $OUTPUT
