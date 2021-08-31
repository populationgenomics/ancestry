#!/usr/bin/env bash

set -ex

gsutil -m mv "gs://cpg-tob-wgs-main-upload/onek1k_scrnaseq_grch38/" $OUTPUT
