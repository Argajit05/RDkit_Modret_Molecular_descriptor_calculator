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

# --- Core Functions ---

def smiles_to_mol(smiles: str):
    """
    Converts a SMILES string to an RDKit Mol object with error handling.

    Args:
        smiles: The SMILES string to convert.

    Returns:
        An RDKit Mol object, or None if conversion or sanitization fails.
    """
    if not isinstance(smiles, str) or "not found" in smiles.lower() or not smiles.strip():
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            Chem.SanitizeMol(mol)
            return mol
        else:
            # MolFromSmiles itself returned None (invalid SMILES)
            logging.debug(f"RDKit returned None for SMILES: {smiles[:50]}...") # Log only beginning
            return None
    except Exception as e:
        logging.warning(f"Error processing SMILES '{smiles[:50]}...': {e}")
        return None

def calculate_rdkit_descriptors(mol):
    """
    Calculates a predefined set of RDKit descriptors for a molecule.

    Args:
        mol: An RDKit Mol object.

    Returns:
        A dictionary mapping descriptor names to their calculated values (or None on error).
    """
    if mol is None:
        return {name: None for name in RDKIT_DESCRIPTOR_FUNCTIONS}

    results = {}
    for name, func in RDKIT_DESCRIPTOR_FUNCTIONS.items():
        try:
            results[name] = func(mol)
        except Exception as e:
            # This catches errors during the calculation for a specific descriptor
            logging.warning(f"Error calculating RDKit descriptor '{name}': {e}")
            results[name] = None
    return results

def calculate_mordred_descriptors(mols_list: list, nproc: int = 1):
    """
    Calculates Mordred descriptors for a list of RDKit Mol objects.

    Args:
        mols_list: A list of valid RDKit Mol objects.
        nproc: Number of processes to use for calculation (default: 1).

    Returns:
        A pandas DataFrame containing Mordred descriptors, indexed like the input list.
        Returns an empty DataFrame if input list is empty or calculation fails globally.
    """
    if not mols_list:
        return pd.DataFrame()

    logging.info(f"Initializing Mordred calculator for {len(descriptors.all)} descriptors...")
    # ignore_3D=True is crucial when working from SMILES without generating 3D coordinates
    mordred_calc = Calculator(descriptors, ignore_3D=True)

    logging.info(f"Calculating Mordred descriptors using {nproc} process(es)...")
    start_time = time.time()

    try:
        # Use mordred's pandas method for efficient bulk calculation
        mordred_results = mordred_calc.pandas(mols_list, nproc=nproc)

        # Convert potential calculation errors (returned as objects) to NaN
        # Ensure numeric types where applicable
        mordred_results = mordred_results.apply(pd.to_numeric, errors='coerce')

        end_time = time.time()
        logging.info(f"Mordred calculation finished in {end_time - start_time:.2f} seconds.")

        # Check for columns that are entirely NaN (might indicate consistent errors)
        all_nan_cols = mordred_results.columns[mordred_results.isna().all()].tolist()
        if all_nan_cols:
             logging.warning(f"{len(all_nan_cols)} Mordred columns have only NaN values and will be kept.")
             # Optional: Drop them if preferred
             # mordred_results = mordred_results.dropna(axis=1, how='all')
             # logging.info("Dropped all-NaN Mordred columns.")

        return mordred_results

    except Exception as e:
        logging.error(f"An error occurred during bulk Mordred calculation: {e}", exc_info=True)
        return pd.DataFrame() # Return empty DataFrame on major failure


# --- Main Execution ---

def main(input_path, output_path, smiles_col='Smiles', n_proc=1):
    """
    Main function to load data, calculate descriptors, and save results.
    """
    logging.info(f"Starting descriptor calculation process.")
    logging.info(f"Input file: {input_path}")
    logging.info(f"Output file: {output_path}")
    logging.info(f"SMILES column: {smiles_col}")
    logging.info(f"Mordred processes: {n_proc}")

    # 1. Load Data
    try:
        df = pd.read_csv(input_path)
        logging.info(f"Loaded dataframe with {len(df)} rows from {input_path}.")
        if smiles_col not in df.columns:
            logging.error(f"SMILES column '{smiles_col}' not found in the input file.")
            sys.exit(1) # Exit if SMILES column is missing
    except FileNotFoundError:
        logging.error(f"Error: Input file not found at '{input_path}'")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error reading CSV file '{input_path}': {e}", exc_info=True)
        sys.exit(1)

    # Keep original index for later merging
    df_original_index = df.index

    # 2. Preprocessing and Molecule Generation
    logging.info("Converting SMILES to RDKit Mol objects...")
    df['mol'] = df[smiles_col].apply(smiles_to_mol)

    original_count = len(df)
    df.dropna(subset=['mol'], inplace=True) # Remove rows where mol generation failed
    valid_count = len(df)
    invalid_count = original_count - valid_count

    if invalid_count > 0:
        logging.warning(f"Removed {invalid_count} rows due to invalid/unparsable SMILES.")
    logging.info(f"Processing {valid_count} valid molecules.")

    if valid_count == 0:
        logging.error("No valid molecules found after SMILES parsing. Exiting.")
        sys.exit(1)

    # 3. RDKit Descriptor Calculation
    logging.info("Calculating RDKit descriptors...")
    start_time = time.time()
    rdkit_desc_list = df['mol'].apply(calculate_rdkit_descriptors).tolist()
    df_rdkit = pd.DataFrame(rdkit_desc_list, index=df.index) # Use valid molecule index
    end_time = time.time()
    logging.info(f"RDKit descriptors calculated in {end_time - start_time:.2f} seconds.")

    # 4. Mordred Descriptor Calculation
    valid_mols_list = df['mol'].tolist() # Get list of valid mol objects
    df_mordred = calculate_mordred_descriptors(valid_mols_list, nproc=n_proc)

    # Ensure Mordred results align with the main DataFrame's valid rows index
    if not df_mordred.empty:
        df_mordred.index = df.index

    # 5. Combine Results
    logging.info("Combining original data with calculated descriptors...")
    # Select original columns (excluding 'mol') and merge based on the index of valid molecules
    df_final = pd.concat([df.drop('mol', axis=1), df_rdkit, df_mordred], axis=1)

    # Reindex to include original rows that might have been dropped due to invalid SMILES, filling descriptors with NaN
    # This ensures the output has the same number of rows as the input, marking invalid SMILES rows clearly.
    # df_final = df_final.reindex(df_original_index) # Uncomment this if you want output rows = input rows
    # logging.info(f"Reindexed final DataFrame to match original {original_count} rows.")


    logging.info(f"Final DataFrame shape: {df_final.shape}")

    # 6. Save Results
    try:
        df_final.to_csv(output_path, index=False)
        logging.info(f"Successfully saved results to {output_path}")
    except Exception as e:
        logging.error(f"Error saving results to '{output_path}': {e}", exc_info=True)
        sys.exit(1)

    logging.info("Descriptor calculation process completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate RDKit and Mordred molecular descriptors from a SMILES file.")
    parser.add_argument("-i", "--input",
                        required=True,
                        help="Path to the input CSV file containing SMILES.")
    parser.add_argument("-o", "--output",
                        required=True,
                        help="Path to save the output CSV file with descriptors.")
    parser.add_argument("-s", "--smiles_column",
                        default="Smiles",
                        help="Name of the column containing SMILES strings (default: Smiles).")
    parser.add_argument("-n", "--nproc",
                        type=int, default=1,
                        help="Number of processes to use for Mordred calculation (default: 1). Use with caution.")

    args = parser.parse_args()

    main(args.input, args.output, args.smiles_column, args.nproc)
