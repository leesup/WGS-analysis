#!/bin/bash

project_name="$1"
project_path="$2"
ref_path="$3"

# Normalize VCFs before Annotation
# Split Multiallelic Sites into Biallelic Records
bcftools norm -m-both "../data/${project_name}.vcf" -o "../data/${project_name}_biallelic.vcf"

# Left-Align and Normalize
bcftools norm -f "${ref_path}/Homo_sapiens_assembly38.fasta" "../data/${project_name}_biallelic.vcf" -o "../data/${project_name}_normalized.vcf"

