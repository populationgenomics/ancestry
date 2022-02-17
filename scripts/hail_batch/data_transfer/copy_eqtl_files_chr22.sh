#!/usr/bin/env bash

set -ex

gsutil cp gs://cpg-tob-wgs-main/scrna-seq/grch38_association_files/genotype_files/tob_genotype_chr22.tsv gs://cpg-tob-wgs-test/scrna_seq/grch38_association_files/genotype_files/
gsutil cp gs://cpg-tob-wgs-main/scrna-seq/grch38_association_files/snp_location_files/snpsloc_chr22.tsv gs://cpg-tob-wgs-test/scrna_seq/grch38_association_files/snp_location_files/
gsutil cp gs://cpg-tob-wgs-main/scrna-seq/grch38_association_files/gene_location_files/GRCh38_geneloc_chr22.tsv gs://cpg-tob-wgs-test/scrna_seq/grch38_association_files/gene_location_files/
