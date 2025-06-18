#!/usr/bin/env python3

import argparse
import json
import os
import time
import pandas as pd
import requests


def parse_arguments():
    """
    Parse command-line arguments for protein data downloader.
    """
    parser = argparse.ArgumentParser(
        description='Download protein data from UniProt'
    )
    parser.add_argument(
        '--input',
        default=('data/1_accessions_input.csv'),
        help='Input CSV file with protein accessions'
    )
    parser.add_argument(
        '--output',
        default=('data/2_protein_dataset.json'),
        help='Output JSON file'
    )
    return parser.parse_args()


def extract_protein_info(data):
    """
    Extract required protein information from UniProt JSON response.

    Args:
        data (dict): JSON data from UniProt API response

    Returns:
        dict: Structured protein information dictionary
    """
    result = {}

    # Extract protein name from description
    protein_desc = data.get("proteinDescription", {})
    recommended = protein_desc.get("recommendedName", {})
    result["protein_name"] = (recommended.get("fullName", {})
                              .get("value", "Unknown"))

    # Primary and secondary accessions
    result["accessions"] = {
        "primary": data.get("primaryAccession", ""),
        "secondary": data.get("secondaryAccessions", [])
    }

    # Gene name extraction
    genes = data.get("genes", [])
    gene_name = "Unknown"
    if genes and "geneName" in genes[0]:
        gene_name = genes[0]["geneName"].get("value", "Unknown")
    result["gene"] = gene_name

    # Entry type/status
    entry_type = data.get("entryType", "")
    result["status"] = entry_type

    # Organism information
    organism = data.get("organism", {})
    result["organism"] = organism.get("scientificName", "Unknown")

    # Natural variants extraction
    variants = []
    for feature in data.get("features", []):
        if feature.get("type") == "Natural variant":
            location = feature.get("location", {})
            start = location.get("start", {}).get("value", "")
            variants.append({
                "position": str(start) if start else "Unknown",
                "description": feature.get("description", "")
            })
    result["variant_ids"] = variants

    # Gene Ontology terms
    go_terms = []
    for ref in data.get("uniProtKBCrossReferences", []):
        if ref.get("database") == "GO":
            go_term = ""
            evidence = ""
            for prop in ref.get("properties", []):
                if prop.get("key") == "GoTerm":
                    go_term = prop.get("value", "")
                elif prop.get("key") == "GoEvidenceType":
                    evidence = prop.get("value", "")

            go_terms.append({
                "go_id": ref.get("id", ""),
                "term": go_term,
                "evidence": evidence
            })
    result["gene_ontology"] = go_terms

    # Disease associations
    diseases = []
    for comment in data.get("comments", []):
        if comment.get("commentType") == "DISEASE":
            disease = comment.get("disease", {})
            diseases.append({
                "diseaseId": disease.get("diseaseId", ""),
                "description": disease.get("description", ""),
                "UniProt classification": disease.get("diseaseAccession", "")
            })
    result["diseases"] = diseases

    # Protein families and domains
    families = []
    domains = []

    # Extract from cross-references (Pfam, Gene3D, SUPFAM)
    for ref in data.get("uniProtKBCrossReferences", []):
        db = ref.get("database", "")
        if db in ["Pfam", "Gene3D", "SUPFAM"]:
            entry_name = ""
            for prop in ref.get("properties", []):
                if prop.get("key") == "EntryName":
                    entry_name = prop.get("value", "")
                    break

            info = {
                "id": ref.get("id", ""),
                "name": entry_name,
                "database": db
            }

            # Classify by database type
            if db in ["Pfam", "SUPFAM"]:
                families.append(info)
            elif db == "Gene3D":
                domains.append(info)

    # Extract domain features
    for feature in data.get("features", []):
        if feature.get("type") == "Domain":
            location = feature.get("location", {})
            domains.append({
                "type": feature.get("type", ""),
                "description": feature.get("description", ""),
                "start": location.get("start", {}).get("value", ""),
                "end": location.get("end", {}).get("value", "")
            })

    result["families_domains"] = {
        "families": families,
        "domains": domains
    }

    # Sequence information and isoforms
    sequence_info = {
        "canonical": data.get("sequence", {}).get("value", ""),
        "isoforms": []
    }

    # Extract alternative isoforms
    for comment in data.get("comments", []):
        if comment.get("commentType") == "ALTERNATIVE PRODUCTS":
            for isoform in comment.get("isoforms", []):
                name_obj = isoform.get("name", {})
                iso_name = name_obj.get("value", "")
                iso_ids = isoform.get("isoformIds", [])

                sequence_info["isoforms"].append({
                    "name": iso_name,
                    "id": iso_ids[0] if iso_ids else "",
                    "status": isoform.get("isoformSequenceStatus", "")
                })

    result["sequence_isoforms"] = sequence_info

    # Functional regions
    regions = []
    for feature in data.get("features", []):
        if feature.get("type") == "Region":
            location = feature.get("location", {})
            regions.append({
                "description": feature.get("description", ""),
                "start": location.get("start", {}),
                "end": location.get("end", {})
            })
    result["regions"] = regions

    return result


def main():
    """
    Main function to execute protein data download pipeline.
    """
    args = parse_arguments()

    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # Read protein accessions from input CSV
    print(f"Reading accessions from: {args.input}")
    df = pd.read_csv(args.input)
    accessions = df['Accession'].str.strip().tolist()
    print(f"Processing {len(accessions)} proteins...")

    # Download protein data from UniProt
    dataset = {}
    for i, acc in enumerate(accessions, 1):
        print(f"Fetching {i}/{len(accessions)}: {acc}")

        try:
            # Make API request to UniProt
            response = requests.get(
                f"https://rest.uniprot.org/uniprotkb/{acc}.json",
                timeout=30
            )

            if response.status_code == 200:
                data = response.json()
                dataset[acc] = extract_protein_info(data)
                print(f"Success: {acc}")
            else:
                print(f"Failed: {acc} (HTTP {response.status_code})")

        except requests.exceptions.RequestException as e:
            print(f"Request failed for {acc}: {e}")

        # Rate limiting to avoid overwhelming the API
        time.sleep(1)

    # Save results to JSON file
    print(f"Saving results to: {args.output}")
    with open(args.output, 'w') as f:
        json.dump(dataset, f, indent=2)

    print(f"\nCompleted! Saved {len(dataset)} proteins to {args.output}")


if __name__ == "__main__":
    main()