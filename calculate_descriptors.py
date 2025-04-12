# 1. Install RDKit, Mordred, and Pandas
# ----------------------------------------
print("Installing RDKit, Mordred, and Pandas...")
!pip install rdkit mordred pandas -q
print("Installation complete.")

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen, MolSurf
from mordred import Calculator, descriptors
import numpy as np
import argparse
import logging
import sys
import time
import warnings

# --- Configuration ---
# Suppress specific warnings if desired (optional)
warnings.filterwarnings("ignore", category=UserWarning, module='mordred')
# from rdkit import rdBase
# rdBase.DisableLog('rdApp.warning') # More aggressive RDKit warning suppression

LOG_FORMAT = '%(asctime)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=LOG_FORMAT, stream=sys.stdout)

# Define which RDKit descriptors to calculate
RDKIT_DESCRIPTOR_FUNCTIONS = {
    'MolWt': Descriptors.MolWt,
    'LogP': Crippen.MolLogP,
    'TPSA': MolSurf.TPSA,
    'NumHDonors': Lipinski.NumHDonors,
    'NumHAcceptors': Lipinski.NumHAcceptors,
    'NumRotatableBonds': Descriptors.NumRotatableBonds,
    'NumRings': Descriptors.RingCount,
    'MolFormula': Descriptors.rdMolDescriptors.CalcMolFormula
}

# 2. Get them to work

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen, MolSurf
from mordred import Calculator, descriptors
import numpy as np
import argparse
import logging
import sys
import time
import warnings

# --- Configuration ---
# Suppress specific warnings if desired (optional)
warnings.filterwarnings("ignore", category=UserWarning, module='mordred')
# from rdkit import rdBase
# rdBase.DisableLog('rdApp.warning') # More aggressive RDKit warning suppression

LOG_FORMAT = '%(asctime)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=LOG_FORMAT, stream=sys.stdout)

# Define which RDKit descriptors to calculate
RDKIT_DESCRIPTOR_FUNCTIONS = {
    'MolWt': Descriptors.MolWt,
    'LogP': Crippen.MolLogP,
    'TPSA': MolSurf.TPSA,
    'NumHDonors': Lipinski.NumHDonors,
    'NumHAcceptors': Lipinski.NumHAcceptors,
    'NumRotatableBonds': Descriptors.NumRotatableBonds,
    'NumRings': Descriptors.RingCount,
    'MolFormula': Descriptors.rdMolDescriptors.CalcMolFormula
}

# 3. Upload the SMILES directly makes it easier! in Goggle Collab.
# ----------------------------------------
csv_content = """Smiles
COCCCCC(=NOCCN)C1=CC=C(C=C1)C(F)(F)F
CCCC(C)(COC(=O)N)COC(=O)N
CN(C)CCC=C1C2=CC=CC=C2CCC3=CC=CC=C31
[Li+].[Li+].[Li+].C(C(=O)[O-])C(CC(=O)[O-])(C(=O)[O-])O
CN1CCN(CC1)C2=NC3=C(C=CC(=C3)Cl)NC4=CC=CC=C42
C1=CC=C(C(=C1)C2=NC(C(=O)NC3=C2C=C(C=C3)Cl)O)Cl
C1C2=NN=CN2C3=C(C=C(C=C3)Cl)C(=N1)C4=CC=CC=C4
CC(C)CC(CC(=O)O)CN
CN=C1CN(C(=C2C=C(C=CC2=N1)Cl)C3=CC=CC=C3)O
CC1=NN=C2N1C3=C(C=C(C=C3)Cl)C(=NC2)C4=CC=CC=C4
CN1CCN(CC1)C2=NC3=CC=CC=C3OC4=C2C=C(C=C4)Cl
CCN(CC)CCN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3F
CNCCCC1C2=CC=CC=C2C=CC3=CC=CC=C13
C1CN(CCN1CCCCC2=CNC3=C2C=C(C=C3)C#N)C4=CC5=C(C=C4)OC(=C5)C(=O)N
CCN(CC)C(=O)C1(CC1CN)C2=CC=CC=C2
CCCCCCCCCCCC(=O)OCN1C(=O)CCC2=C1C=C(C=C2)OCCCCN3CCN(CC3)C4=C(C(=CC=C4)Cl)Cl
CCCCCCC(C)(C)C1=CC(=C2C3CC(=O)CCC3C(OC2=C1)(C)C)O
C1CCC(C(C1)CN2CCN(CC2)C3=NSC4=CC=CC=C43)CN5C(=O)C6C7CCC(C7)C6C5=O
CN1C(=O)CC(=O)N(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3
CN(C)CCCC1(C2=C(CO1)C=C(C=C2)C#N)C3=CC=C(C=C3)F
C1CNCC(C1C2=CC=C(C=C2)F)COC3=CC4=C(C=C3)OCO4
C1CN(CCN1CCC2=C(C=C3C(=C2)CC(=O)N3)Cl)C4=NSC5=CC=CC=C54
CC1=CC2=C(S1)NC3=CC=CC=C3N=C2N4CCN(CC4)C
C1CN(CCC1(C2=CC=C(C=C2)Cl)O)CCCC(=O)C3=CC=C(C=C3)F
C1CN(CCN1)C2=NC3=CC=CC=C3OC4=C2C=C(C=C4)Cl
C1CN(CCN1CCCN2C(=O)N3C=CC=CC3=N2)C4=CC(=CC=C4)Cl
CCOC1=CC=CC=C1OCC2CNCCO2
CC(CN1C2=CC=CC=C2CCC3=CC=CC=C31)CN(C)C
CNCCCN1C2=CC=CC=C2CCC3=CC=CC=C31
CCCC(CCC)C(=O)O
CNCCC=C1C2=CC=CC=C2CCC3=CC=CC=C31
CN(C)CCCN1C2=CC=CC=C2CCC3=CC=CC=C31
CC1C(OCCN1C)C2=CC=CC=C2
CC1=NN=C2N1C3=C(C=C(C=C3)Cl)C(=NC2)C4=CC=CC=C4Cl
CCCC(C)(COC(=O)N)COC(=O)NC(C)C
CC1(CCC2C(C1)CCC3C2CCC4(C3CCC4C(=O)CN5C=C(C=N5)C#N)C)O
C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)N
CC1(CC(=O)N(C(=O)C1)CCCCN2CCN(CC2)C3=NC=CC=N3)C
C1=CC=C(C=C1)C2=NC(C(=O)NC3=C2C=C(C=C3)Cl)O
C1CCC(CC1)(CC(=O)O)CN
CCC1=C(NC2=C1C(=O)C(CC2)CN3CCOCC3)C
[Li+].[Li+].C(=O)([O-])[O-]
C1CN(CC=C1N2C3=CC=CC=C3NC2=O)CCCC(=O)C4=CC=C(C=C4)F
C1=CC(=C(C(=C1)Cl)Cl)C2=C(N=C(N=N2)N)N
CCN1CCCC1CNC(=O)C2=CC(=C(C=C2OC)N)S(=O)(=O)CC
CC1=NC=C2N1C3=C(C=C(C=C3)Cl)C(=NC2)C4=CC=CC=C4F
CNC1(CCCCC1=O)C2=CC=CC=C2Cl
"""

# Save the string content to a file in Colab's temporary storage
csv_filename = "ALLDrugswithSMILES.csv"
with open(csv_filename, "w") as f:
    f.write(csv_content)

# Now load the DataFrame
try:
    df = pd.read_csv(csv_filename)
    print(f"\nLoaded dataframe with {len(df)} rows from {csv_filename}.")
    print("DataFrame head:")
    print(df.head())
except FileNotFoundError:
    print(f"Error: File '{csv_filename}' not found. Make sure it was created or uploaded correctly.")
    # Stop execution if file loading fails
    raise SystemExit("Stopping execution due to missing file.")
except Exception as e:
    print(f"An error occurred while reading the CSV: {e}")
    raise SystemExit("Stopping execution due to CSV reading error.")

# 4. Converting SMILES to RDKit Mol objects
# ----------------------------------------
def smiles_to_mol(smiles):
    try:
        # Handle potential non-string inputs or common placeholders
        if not isinstance(smiles, str) or "not found" in smiles.lower() or smiles.strip() == "":
             return None
        mol = Chem.MolFromSmiles(smiles)
        # Basic sanitization check (optional but good practice)
        if mol:
             Chem.SanitizeMol(mol)
        return mol # Returns None if Chem.MolFromSmiles fails
    except Exception as e:
        # print(f"Error processing SMILES '{smiles}': {e}") # Uncomment for debugging
        return None

# Apply the conversion
print("\nConverting SMILES to RDKit Mol objects...")
df['mol'] = df['Smiles'].apply(smiles_to_mol)

# Check for and remove rows with invalid SMILES
original_count = len(df)
df.dropna(subset=['mol'], inplace=True)
valid_count = len(df)
invalid_count = original_count - valid_count

if invalid_count > 0:
    print(f"Warning: Removed {invalid_count} rows due to invalid or unparsable SMILES.")
print(f"Processing {valid_count} valid molecules.")

# 5. RDKit Descriptor Calculation
# ----------------------------------------
print("\nCalculating RDKit descriptors...")

def calculate_rdkit_descriptors(mol):
    """Calculates a selected set of RDKit descriptors."""
    if mol is None:
        # Return a dictionary with None values for consistency
        return {'MolWt': None, 'LogP': None, 'TPSA': None, 'NumHDonors': None,
                'NumHAcceptors': None, 'NumRotatableBonds': None, 'NumRings': None,
                'MolFormula': None}
    try:
        return {
            'MolWt': Descriptors.MolWt(mol),
            'LogP': Crippen.MolLogP(mol), # Using Crippen's LogP calculation
            'TPSA': MolSurf.TPSA(mol),    # Using MolSurf's TPSA calculation
            'NumHDonors': Lipinski.NumHDonors(mol),
            'NumHAcceptors': Lipinski.NumHAcceptors(mol),
            'NumRotatableBonds': Descriptors.NumRotatableBonds(mol),
            'NumRings': Descriptors.RingCount(mol),
            'MolFormula': Descriptors.rdMolDescriptors.CalcMolFormula(mol) # Example: Calculate Mol Formula
        }
    except Exception as e:
        print(f"Error calculating RDKit descriptors for a molecule: {e}")
        # Return None dict in case of calculation error for a valid mol object
        return {'MolWt': None, 'LogP': None, 'TPSA': None, 'NumHDonors': None,
                'NumHAcceptors': None, 'NumRotatableBonds': None, 'NumRings': None,
                'MolFormula': None}

# Apply the function and create a new DataFrame
rdkit_desc_list = df['mol'].apply(calculate_rdkit_descriptors).tolist()
df_rdkit = pd.DataFrame(rdkit_desc_list, index=df.index)

print("RDKit descriptors calculated.")
print("RDKit DataFrame head:")
print(df_rdkit.head())

# 6. Mordred Descriptor Calculation
# ----------------------------------------
print("\nCalculating Mordred descriptors (this may take some time)...")

# Initialize Mordred calculator
# ignore_3D=True is crucial when working from SMILES without generating 3D coordinates
mordred_calc = Calculator(descriptors, ignore_3D=True)

# Mordred can calculate descriptors in bulk, which is faster than row-by-row apply
# Filter out any None molecules just in case (should be handled by dropna already)
valid_mols = df['mol'].dropna().tolist()

if valid_mols:
    # Calculate descriptors for all valid molecules
    # The result is a pandas Series/DataFrame indexed like the input list
    try:
        mordred_results = mordred_calc.pandas(valid_mols, nproc=1) # nproc=1 avoids potential multiprocessing issues in Colab

        # Mordred might return results with mixed types (e.g., errors as objects)
        # Convert errors to NaN and ensure numeric types where possible
        mordred_results = mordred_results.apply(pd.to_numeric, errors='coerce')

        # Set the index of mordred results to match the original DataFrame's valid rows
        mordred_results.index = df.index # Crucial for correct alignment

        # Fill missing values (e.g., for molecules where calculation failed) with NaN
        df_mordred = mordred_results.reindex(df.index)

        print(f"Mordred descriptors calculated ({df_mordred.shape[1]} descriptors).")
        print("Mordred DataFrame head:")
        print(df_mordred.head())

        # Optional: Check for columns that are entirely NaN (might indicate consistent errors)
        all_nan_cols = df_mordred.columns[df_mordred.isna().all()].tolist()
        if all_nan_cols:
             print(f"\nWarning: {len(all_nan_cols)} Mordred columns have only NaN values.")
             # print("All NaN columns:", all_nan_cols) # Uncomment to see the names
             # df_mordred = df_mordred.dropna(axis=1, how='all') # Optionally drop them
             # print("Dropped all-NaN columns.")

    except Exception as e:
        print(f"An error occurred during bulk Mordred calculation: {e}")
        # Create an empty DataFrame as a fallback
        df_mordred = pd.DataFrame(index=df.index)
else:
    print("No valid molecules found for Mordred calculation.")
    df_mordred = pd.DataFrame(index=df.index)

# 7. Combine Results
# ----------------------------------------
print("\nCombining original data with calculated descriptors...")

# Drop the temporary 'mol' column from the original df before concatenating
df_combined = pd.concat([df.drop('mol', axis=1), df_rdkit, df_mordred], axis=1)

print("Final DataFrame shape:", df_combined.shape)
print("Final DataFrame head:")
print(df_combined.head())
print("\nFinal DataFrame info:")
df_combined.info()

# 8. Save to CSV
# ----------------------------------------
output_filename = "ALLDrugswithSMILES_and_Descriptors.csv"
df_combined.to_csv(output_filename, index=False)
print(f"\nResults saved to {output_filename}")

# 9. Offer download link in Colab
from google.colab import files
files.download(output_filename)
