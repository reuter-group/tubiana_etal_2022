#!/usr/bin/env python3
"""
PDB Processing and Protrusions Calculator

This script processes PDB files and computes protrusions using the PePr2DS pipeline.
It corresponds to steps 2 and 3 from the original notebook workflow.

Usage:
    python process_pdbs_and_protrusions.py <input_folder> <domain_name> [--output_csv <path>]

Example:
    python process_pdbs_and_protrusions.py ./pdb_files START --output_csv ./output.csv
"""

import argparse
import glob
import os
import tempfile
import sys
from pathlib import Path
import pandas as pd
from urllib import request
import gzip
from io import BytesIO
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from tqdm import tqdm
from biopandas.pdb import PandasPdb

def download_mapping_file(tmpdir):
    """Download the PDB chain uniprot mapping file"""
    url = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_uniprot.csv.gz"
    out_file_path = url.split("/")[-1][:-3]
    
    print(f'Downloading pdb_chain_uniprot.csv Database from: {url}')
    
    try:
        response = request.urlopen(url)
        compressed_file = BytesIO(response.read())
        decompressed_file = gzip.GzipFile(fileobj=compressed_file, mode='rb')

        # Extract database
        output_path = os.path.join(tmpdir, out_file_path)
        with open(output_path, 'w') as outfile:
            outfile.write(decompressed_file.read().decode('utf-8'))

        return output_path
    except Exception as e:
        print(f"Error downloading mapping file: {e}")
        return None

def process_pdbs_and_protrusions(input_folder, domain_name, output_csv=None, compute_extra_features=False):
    """
    Process PDB files and compute protrusions
    
    Args:
        input_folder (str): Folder containing PDB files
        domain_name (str): Domain name for the dataset (max 4 letters recommended)
        output_csv (str): Output CSV file path (optional)
        compute_extra_features (bool): Whether to compute extra features (dssp, SASA etc...) (default: False)

    Returns:
        pandas.DataFrame: Dataset with protrusions computed
    """
    
    # Validate inputs
    if not os.path.exists(input_folder):
        raise ValueError(f"Input folder does not exist: {input_folder}")
    
    # Find PDB files
    pdb_pattern = os.path.join(input_folder, "*.pdb")
    pdblist = glob.glob(pdb_pattern)
    
    if not pdblist:
        raise ValueError(f"No PDB files found in {input_folder}")
    
    print(f"Found {len(pdblist)} PDB files in {input_folder}")
    
    # Create temporary directory
    tmpdir = tempfile.mkdtemp()
    
    try:
        # Import PePr2DS modules
        print("Initializing PePr2DS builder...")
        import importlib
        import pepr2ds.builder.Builder as builderEngine
        importlib.reload(builderEngine)
        
        # Setup configuration dictionary
        SETUP = {
            "DOMAIN_PROSITE": None,
            "PROSITE_DOMAIN": None,
            "DOMAIN_CATH": None,
            "CATH_DOMAIN": None,
            "SUPERFAMILY": None
        }
        
        # Initialize builder
        builder = builderEngine.Builder(
            SETUP, 
            recalculate=True, 
            update=False, 
            notebook=False,  # Changed to False since we're not in a notebook 
            core=1
        )
        
        # Configure builder
        builder.structure.CATHFOLDER = tmpdir  # Use tmpdir instead of output_folder
        builder.structure.ALFAFOLDFOLDER = ''
        
      
        
        if compute_extra_features:

            # Download mapping file
            print("Downloading PDB chain mapping file...")
            mappingFile = download_mapping_file(tmpdir)
            if mappingFile is None:
                raise RuntimeError("Failed to download mapping file")
            
            # Step 2: Process PDB files
            print("Step 2: Processing PDB files...")
            print("This may take a while depending on the number and size of PDB files...")

            DATASET = builder.structure.process_pdb_list(
                pdblist=pdblist, 
                datatype='custom',
                domname=domain_name, 
                mappingFile=mappingFile
            )
        
            if DATASET is None or DATASET.empty:
                raise RuntimeError("Failed to process PDB files - no data generated")
        
        else: 
            #Read just all the pdb from pdblist with biopandas
            datasetlist = []
            for pdb in pdblist:
                pdbdf = PandasPdb().read_pdb(pdb).df["ATOM"]
                pdbdf['cathpdb'] = pdb.split("/")[-1].split(".")[0]
                pdbdf["domain"] = domain_name
                datasetlist.append(pdbdf)

            DATASET = pd.concat(datasetlist, ignore_index=True)
        

        print(f"Processed {len(DATASET)} records from PDB files")
        
        # Step 3: Compute protrusions
        print("Step 3: Computing protrusions...")
        print("This step may also take some time...")
        
        


        DATASET = builder.structure.add_protrusions(DATASET)

        
        if DATASET is None or DATASET.empty:
            raise RuntimeError("Failed to compute protrusions")
        
        print(f"Computed protrusions for {len(DATASET)} records")
        
        # Remove duplicates
        print("Removing duplicates...")
        initial_count = len(DATASET)
        DATASET = DATASET.drop_duplicates(
            subset=['atom_number', 'atom_name', 'residue_name', 'residue_number', 'cathpdb', 'chain_id']
        )
        final_count = len(DATASET)
        print(f"Removed {initial_count - final_count} duplicate records")
        # Save output CSV if specified
        #Keep only CB atoms (where protrusions are defined) and CA of GLY amino acids (that don't have CB)
        DATASET = DATASET[
            (DATASET['atom_name'] == 'CB') |
            ((DATASET['residue_name'] == 'GLY') & (DATASET['atom_name'] == 'CA'))
        ]
        if output_csv:
            print(f"Saving dataset to {output_csv}")
            DATASET.to_csv(output_csv, index=False)
        
        print("Processing completed successfully!")
        return DATASET
        
    except ImportError as e:
        print(f"Error importing PePr2DS modules: {e}")
        print("Make sure PePr2DS is properly installed and accessible")
        raise
    except Exception as e:
        print(f"Error during processing: {e}")
        raise
    finally:
        # Cleanup temporary directory
        import shutil
        try:
            shutil.rmtree(tmpdir)
        except:
            pass

def main():
    """Main function for command line interface"""
    parser = argparse.ArgumentParser(
        description="Process PDB files and compute protrusions using PePr2DS",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python process_pdbs_and_protrusions.py ./pdb_files START --output_csv ./output.csv
  python process_pdbs_and_protrusions.py ./pdb_files START -v --output_csv ./results.csv
        """
    )
    
    parser.add_argument(
        'input_folder',
        help='Folder containing PDB files to process'
    )
    
    parser.add_argument(
        '--domain_name', "-n",
        help='Domain name for the dataset (default is input folder)',
        default=None,
    )

    parser.add_argument(
        '--output_csv',
        help='Output CSV file path (required if you want to save the results, default is domain name)',
        default=None,
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output'
    )

    parser.add_argument(
        '--compute_extra_features',
        action='store_true',
        help='Compute extra features (default is False)'
    )

    args = parser.parse_args()

    #If a domain name is given, takes the domain name, otherwise take the name of the input folder
    args.domain_name = args.domain_name if args.domain_name else args.input_folder
    args.output_csv = args.output_csv if args.output_csv else f"{args.domain_name}.csv"
    
    # Suppress BioPython warnings unless verbose
    if not args.verbose:
        warnings.filterwarnings('ignore', category=PDBConstructionWarning)
    
    try:
        # Validate domain name
        if len(args.domain_name) > 4:
            print(f"Warning: Domain name '{args.domain_name}' is longer than 4 characters")
        
        # Process PDBs and compute protrusions
        dataset = process_pdbs_and_protrusions(
            input_folder=args.input_folder,
            domain_name=args.domain_name,
            output_csv=args.output_csv,
            compute_extra_features=args.compute_extra_features
        )
        
        print(f"\nSummary:")
        print(f"- Input folder: {args.input_folder}")
        print(f"- Domain name: {args.domain_name}")
        print(f"- Records processed: {len(dataset)}")
        #print(f"- Columns: {list(dataset.columns)}")
        
        
        print(f"- Output CSV: {args.output_csv}")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
