# PePrMInt Utilities

This directory contains three main utilities for working with PDB structures and generating custom datasets for PePrMInt analysis. These tools provide a complete workflow from downloading structures to computing structural features and generating final datasets.

## Prerequisites

Before using these tools, you need to set up the PePrMInt environment:

1. **Install the PePrMInt environment:**
   ```bash
   conda env create -f environment.yml
   conda activate peprmint
   ```

2. **Install PePr2DS in development mode:**
   ```bash
   cd /home/thibault/projects/peprmint/tubiana_etal_2022
   pip install -e .
   ```

## Tools Overview

### 1. 📥 `download_pdb_from_pfam.py` - PDB Structure Downloader

Downloads all PDB structures containing a specific Pfam domain family. This tool queries the RCSB PDB database and downloads relevant structures for further analysis.

**Usage:**
```bash
python download_pdb_from_pfam.py <PFAM_ID> [options]
```

**Example:**
```bash
# Download all structures containing the START domain (PF01852)
python download_pdb_from_pfam.py PF01852

# Download to a custom folder
python download_pdb_from_pfam.py PF01852 --output START_domain_structures

# Download RING finger domain structures
python download_pdb_from_pfam.py PF00097 --output RING_structures
```

**Features:**
- Automatic querying of RCSB PDB database
- Fallback to UniProt database if RCSB search fails
- Skips already downloaded files
- Rate-limited downloads to be respectful to servers
- Creates organized output directories

### 2. 🔬 `process_pdbs_and_protrusions.py` - Structural Feature Calculator

Computes comprehensive structural features for PDB files including SASA (Solvent Accessible Surface Area), secondary structure, and protrusions. This tool processes the structural data without computing the interfacial binding sites (IBS)

**Usage:**
```bash
python process_pdbs_and_protrusions.py <input_folder> [options]
```

**Examples:**
```bash
# Basic usage - process all PDB files in a folder
python process_pdbs_and_protrusions.py ./START_structures --domain_name START

# Save results to specific CSV file
python process_pdbs_and_protrusions.py ./RING_structures --domain_name RING --output_csv ring_features.csv

# Enable verbose output and compute extra features
python process_pdbs_and_protrusions.py ./structures --domain_name CUSTOM -v --compute_extra_features
```

**Output:**
- CSV file containing all computed features
- Focus on CB atoms (CA for glycine residues)
- Duplicate removal and data cleaning

### 3. 🌐 `generate_custom_csv.ipynb` - Interactive Web Interface

A Jupyter notebook designed to work with Voila, providing an interactive web interface for generating custom PePrMInt datasets. This tool allows users to define interfacial binding sites (IBS) according to the methodology described in Tubiana et al. 2022.

**Installation & Usage:**

1. **Start the web interface:**
   ```bash
   voila generate_custom_csv.ipynb
   ```

2. **Open your web browser** and navigate to the localhost URL displayed in the terminal (typically `http://localhost:8866/`)

**Features:**

- **📁 Folder Selection**: Choose input folder containing PDB files and output directory
- **⚙️ Alignment Options**: 
  - TMalign structural alignment (recommended)
  - Option to use pre-aligned structures
- **🎯 Reference Structure Definition**:
  - Select reference PDB file
  - Define 3 amino acid positions for binding site orientation
  - Real-time 3D visualization for validation
- **📊 IBS Computation**: Automatic calculation of interfacial binding sites based on geometric criteria
- **💾 Export**: Generate final CSV files compatible with PePrMInt analysis


## Citation

If you use these tools in your research, please cite:

> Tubiana T, Sillitoe I, Orengo C, Reuter N (2022) Dissecting peripheral protein-membrane interfaces. PLOS Computational Biology 18(12): e1010346. https://doi.org/10.1371/journal.pcbi.1010346

