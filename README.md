# Molecular Descriptor Calculator

Calculates a comprehensive set of 2D molecular descriptors using RDKit and Mordred from a list of SMILES strings.

## Features

*   Reads SMILES from a CSV file.
*   Handles invalid SMILES strings gracefully.
*   Calculates a standard set of RDKit descriptors (e.g., MolWt, LogP, TPSA, H-bond donors/acceptors, Rotatable Bonds, Rings, Formula).
*   Calculates a wide range of 2D descriptors using the Mordred library.
*   Outputs the original data merged with the calculated descriptors to a new CSV file.
*   Provides logging for progress and warnings.
*   Configurable via command-line arguments.

## Prerequisites

*   Python (>= 3.8 recommended)
*   Pip (Python package installer)

## Installation

1.  **Clone the repository:**
    ```bash
    git clone <your-repository-url>
    cd molecular_descriptor_calculator
    ```

2.  **Create and activate a virtual environment (Recommended):**
    ```bash
    python -m venv venv
    # On Windows:
    # venv\Scripts\activate
    # On macOS/Linux:
    # source venv/bin/activate
    ```

3.  **Install dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

## Usage

Run the script from the command line, providing the input and output file paths.

**Basic Usage:**

```bash
python calculate_descriptors.py -i data/ALLDrugswithSMILES.csv -o data/output_descriptors.csv
```
## Options:
-i or --input: Path to the input CSV file (required).
-o or --output: Path for the output CSV file (required).
-s or --smiles_column: Name of the column containing SMILES (default: Smiles).
-n or --nproc: Number of processes for Mordred calculation (default: 1). Note: Using multiple processes (nproc > 1) might not always speed things up significantly due to overhead and can sometimes cause issues in certain environments.
Example with custom SMILES column:
python calculate_descriptors.py -i data/my_molecules.csv -o data/my_molecules_descriptors.csv -s Canonical_SMILES

## Input File Format
The input file must be a CSV file containing at least one column with SMILES strings. By default, the script looks for a column named Smiles. Use the -s argument if your column has a different name.
```
Example data/ALLDrugswithSMILES.csv:
Smiles,OtherData
COCCCCC(=NOCCN)C1=CC=C(C=C1)C(F)(F)F,Value1
CCCC(C)(COC(=O)N)COC(=O)N,Value2
```
## Output File Format
The output CSV file will contain all columns from the input file, plus new columns for each calculated RDKit and Mordred descriptor. Rows corresponding to invalid input SMILES will have NaN values in the descriptor columns.
```
Example data/output_descriptors.csv:
Smiles,OtherData,MolWt,LogP,TPSA,...,mordred_descriptor_N
COCCCCC(=NOCCN)C1=CC=C(C=C1)C(F)(F)F,Value1,378.36,3.95,58.12,...,12.34
CCCC(C)(COC(=O)N)COC(=O)N,Value2,248.28,0.45,116.58,...,5.67
```` 
