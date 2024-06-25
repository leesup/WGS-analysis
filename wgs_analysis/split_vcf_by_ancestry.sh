#!/bin/bash

project_name="$1"
project_path="$2"

# Split VCF by Ancestry
bcftools view -S "${project_path}/${project_name}_ids.txt" --force-samples "${project_path}/joint-genotyping.vcf.gz" -o "../data/${project_name}.vcf"

