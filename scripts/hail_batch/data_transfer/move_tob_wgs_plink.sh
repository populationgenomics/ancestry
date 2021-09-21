#!/usr/bin/env bash

set -ex

gsutil -m -cp -r "gs://cpg-tob-wgs-main/tob_wgs_plink/v0/" "gs://cpg-tob-wgs-shared-jpowell/mt/plink/"
