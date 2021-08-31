#!/usr/bin/env bash

set -ex

gsutil -m mv "gs://cpg-tob-wgs-main/sc-rnaseq" "gs://cpg-tob-wgs-main/scrna-seq"
