#!/usr/bin/env bash

set -ex

gsutil cp -m -r "gs://cpg-tob-wgs-main/tob_wgs_snp_chip_variant_pca/v2/scores.ht" $OUTPUT
