import requests
import json
from pathlib import Path
import urllib.request
import time
import argparse

def download_pfam_pdbs(pfam_id, output_folder="pfam_pdbs"):
    """
    Download PDB files for structures containing a specific Pfam domain.
    
    Args:
        pfam_id (str): Pfam family identifier (e.g., "PF01852")
        output_folder (str): Folder to save downloaded PDB files
    """
    
    # Create output directory
    output_path = Path(output_folder)
    output_path.mkdir(exist_ok=True)
    
    print(f"Searching for PDB structures containing Pfam family {pfam_id}...")
    
    # Query PDBe API for structures containing the Pfam domain
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name",
                "operator": "exact_match",
                "value": "Eukaryota"
            }
        },
        "request_options": {
            "return_all_hits": True
        },
        "return_type": "entry"
    }
    
    # Alternative query specifically for Pfam
    pfam_query = {
        "query": {
            "type": "terminal",
            "service": "text", 
            "parameters": {
                "attribute": "rcsb_polymer_entity_annotation.annotation_id",
                "operator": "exact_match",
                "value": pfam_id
            }
        },
        "request_options": {
            "return_all_hits": True
        },
        "return_type": "entry"
    }
    
    try:
        # Search for structures with the Pfam domain
        response = requests.post(search_url, json=pfam_query, timeout=30)
        response.raise_for_status()
        
        search_results = response.json()
        
        if "result_set" not in search_results:
            print(f"No structures found for Pfam {pfam_id}")
            return
            
        pdb_ids = [item["identifier"] for item in search_results["result_set"]]
        
        if not pdb_ids:
            print(f"No PDB structures found for Pfam {pfam_id}")
            return

        print(f"Found {len(pdb_ids)} structures. Downloading...")
        downloaded = 0
        failed = 0
        
        for pdb_id in pdb_ids:
            try:
                # Download PDB file
                pdb_url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
                output_file = output_path / f"{pdb_id.upper()}.pdb"
                
                if output_file.exists():
                    print(f"  {pdb_id}: already exists, skipping")
                    continue
                
                print(f"  Downloading {pdb_id}...")
                urllib.request.urlretrieve(pdb_url, output_file)
                downloaded += 1
                
                # Be nice to the server
                time.sleep(0.1)
                
            except Exception as e:
                print(f"  Failed to download {pdb_id}: {e}")
                failed += 1
                
        print(f"\nDownload complete!")
        print(f"  Successfully downloaded: {downloaded} files")
        print(f"  Failed downloads: {failed} files")
        print(f"  Files saved to: {output_path.absolute()}")
        
    except requests.RequestException as e:
        print(f"Error searching for Pfam {pfam_id}: {e}")
        print("Trying alternative approach with UniProt...")
        
        # Fallback: Use UniProt API to find proteins with the Pfam domain
        try:
            uniprot_url = f"https://rest.uniprot.org/uniprotkb/search"
            params = {
                "query": f"family:{pfam_id}",
                "format": "json",
                "size": max_structures
            }
            
            uniprot_response = requests.get(uniprot_url, params=params, timeout=30)
            uniprot_response.raise_for_status()
            
            uniprot_data = uniprot_response.json()
            
            if "results" not in uniprot_data:
                print("No UniProt entries found")
                return
                
            print(f"Found {len(uniprot_data['results'])} UniProt entries")
            
            # Extract PDB IDs from UniProt cross-references
            pdb_ids_from_uniprot = []
            for entry in uniprot_data["results"]:
                if "uniProtKBCrossReferences" in entry:
                    for xref in entry["uniProtKBCrossReferences"]:
                        if xref["database"] == "PDB":
                            pdb_ids_from_uniprot.append(xref["id"])
            
            if not pdb_ids_from_uniprot:
                print("No PDB structures found in UniProt cross-references")
                return
                
            print(f"Found {len(set(pdb_ids_from_uniprot))} unique PDB IDs from UniProt")
            
            # Download the PDB files
            downloaded = 0
            for pdb_id in set(pdb_ids_from_uniprot[:max_structures]):
                try:
                    pdb_url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
                    output_file = output_path / f"{pdb_id.upper()}.pdb"
                    
                    if output_file.exists():
                        continue
                        
                    print(f"  Downloading {pdb_id}...")
                    urllib.request.urlretrieve(pdb_url, output_file)
                    downloaded += 1
                    time.sleep(0.1)
                    
                except Exception as e:
                    print(f"  Failed to download {pdb_id}: {e}")
            
            print(f"Downloaded {downloaded} files via UniProt fallback")
            
        except Exception as e:
            print(f"UniProt fallback also failed: {e}")

def main():
    parser = argparse.ArgumentParser(
        description="Download PDB files for structures containing a specific Pfam domain.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        "pfam_id", 
        help="Pfam family identifier (e.g., PF01852)"
    )
    
    parser.add_argument(
        "-o", "--output", 
        help="Output folder to save downloaded PDB files (default: pfam_id)",
        default=None
    )
    
    args = parser.parse_args()
    
    # Use pfam_id as default output folder if not specified
    output_folder = args.output if args.output else args.pfam_id
    
    download_pfam_pdbs(args.pfam_id, output_folder)

if __name__ == "__main__":
    main()