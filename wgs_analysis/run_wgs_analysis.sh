#!/bin/bash

# Load Applications
module load bcftools

# Define variable
project_name="chile"
project_path="/data/CARD/projects/Paula_WGS"
ref_path="/data/CARD/projects/fulyaak"

# Change directory to the location of run_wgs_analysis.sh
cd "$(dirname "$0")"

# Run fetch_gene_info.py to generate gene_dict and capture output
chmod +x fetch_gene_info.py
gene_dict=$(python3 fetch_gene_info.py)

# Check if gene_dict is populated correctly
if [ -z "$gene_dict" ]; then
    echo "Error: gene_dict is empty or not generated correctly."
    exit 1
fi

echo "Printing gene_dict:"
echo "$gene_dict"

# Split VCF by Ancestry
chmod +x split_vcf_by_ancestry.sh
./split_vcf_by_ancestry.sh "$project_name" "$project_path"

# Normalize VCFs before Annotatio
chmod +x normalize_vcf.sh
./normalize_vcf.sh "$project_name" "$project_path" "$ref_path"

