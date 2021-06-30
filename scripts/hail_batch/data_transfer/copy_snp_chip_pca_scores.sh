#!/usr/bin/env bash

set -ex

gsutil -m cp -r "gs://cpg-tob-wgs-main/tob_wgs_snp_chip_variant_pca/v2/scores.ht" $OUTPUT
