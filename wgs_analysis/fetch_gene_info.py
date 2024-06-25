#!/usr/bin/env python3

# Import Python Packages
import requests

def fetch_gene_info_ensembl(file_gene_names="gene_names.txt", species='human', genome_version='GRCh38'):
    gene_names = []

    # Open gene_names.txt located in the same directory
    with open(file_gene_names, 'r') as file:
        for line in file:
            gene_name = line.strip()
            gene_names.append(gene_name)

    gene_info_dict = {}
    server = "https://rest.ensembl.org"

    for gene_name in gene_names:
        endpoint = f"/lookup/symbol/{species}/{gene_name}"
        headers = {"Content-Type": "application/json"}

        try:
            response = requests.get(server + endpoint, headers=headers, params={"expand": "1"})
            response.raise_for_status()  # Raise HTTPError for bad responses
        except requests.exceptions.RequestException as e:
            print(f"Fetching failed for {gene_name}: {str(e)}")
            continue

        data = response.json()
        gene_info = {
            "gene_name": data.get("display_name", gene_name),
            "chromosome": f"chr{data['seq_region_name']}",
            "start": int(data["start"]),
            "end": int(data["end"]),
            "genome_version": genome_version
        }

        gene_info_dict[gene_name] = gene_info

    return gene_info_dict

if __name__ == "__main__":
    gene_dict = fetch_gene_info_ensembl()
    print(gene_dict)  # Print gene_dict for testing/debugging

