#!/bin/bash

set -ex

# Copy a subset of the data to the `test` bucket.
gsutil -m cp -r gs://cpg-tob-wgs-main/joint_vcf/v1/raw/genomes.mt gs://cpg-tob-wgs-test/joint_vcf/v1/raw/
