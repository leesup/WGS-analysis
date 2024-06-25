#!/bin/bash

# Path of Data Directories
path_data = "/data/CARD/projects/Paula_WGS"
path_ref = "/data/CARD/projects/fulyaak"

# Change directory to the location of run_wgs_analysis.sh
cd "$(dirname "$0")"

# Ensure fetch_gene_info.py has executable permissions
chmod +x fetch_gene_info.py

# Run fetch_gene_info.py to generate gene_dict and capture output
gene_dict=$(python3 fetch_gene_info.py)

# Check if gene_dict is populated correctly
if [ -z "$gene_dict" ]; then
    echo "Error: gene_dict is empty or not generated correctly."
    exit 1
fi

echo "Printing gene_dict:"
echo "$gene_dict"

# Load Applications on Biowulf
module load bcftools

# Split VCF by Ancestry
bcftools view -S "${path_data}/brasil_ids.txt" "${path_data}/brasil.vcf" "${path_data}/joint-genotyping.vcf.gz"
bcftools view -S "${path_data}/chile_ids.txt" "${path_data}/chile.vcf" "${path_data}/joint-genotyping.vcf.gz"

# Normalize VCFs before Annotation
# Split Multiallelic Sites into Biallelic Records
bcftools norm -m-both -o "${path_data}/brasil_biallelic.vcf" "${path_data}/brasil.vcf"
bcftools norm -m-both -o "${path_data}/chile_biallelic.vcf" "${path_data}/chile.vcf"

# Left-Align and Normalize
bcftools norm -f "${path_ref}/Homo_sapiens_assembly38.fasta" -o "${path_data}/brasil_normalized.vcf" "${work_dir}/brasil_biallelic.vcf"
bcftools norm -f "${path_Ref}/Homo_sapiens_assembly38.fasta" -o "${path_data}/chile_normalized.vcf" "${work_dir}/chile_biallelic.vcf"

