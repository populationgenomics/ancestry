#!/usr/bin/env python3

"""
Copies PCA scores from TOB-WGS/SNP-chip data
to test bucket.
"""

import os
import hail as hl

output = os.getenv('OUTPUT')

INPUT_MT = 'gs://cpg-tob-wgs-main/tob_wgs_snp_chip_variant_pca/v2/scores.ht'

mt = hl.read_table(INPUT_MT)
output_path = f'{output}/scores.ht'
mt.head(n = 10000).write(output_path)
