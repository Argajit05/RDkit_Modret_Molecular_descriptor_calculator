
# Troubleshooting Guide

This guide covers common issues encountered when running the `calculate_descriptors.py` script.

## 1. Installation Errors

*   **`ModuleNotFoundError: No module named 'rdkit'` (or `mordred`, `pandas`)**
    *   **Cause:** The required library is not installed in the current Python environment.
    *   **Solution:**
        1.  Ensure you have activated your virtual environment (if you created one): `source venv/bin/activate` (Linux/macOS) or `venv\Scripts\activate` (Windows).
        2.  Run `pip install -r requirements.txt` again to install all dependencies.
        3.  Verify the installation: `pip show rdkit mordred pandas numpy`.

*   **Errors during RDKit installation:**
    *   **Cause:** RDKit installation can sometimes be complex due to its C++ dependencies. Using `pip install rdkit` is generally the most reliable method now. Conda (`conda install -c conda-forge rdkit`) used to be preferred but pip support has improved significantly.
    *   **Solution:** Ensure you have necessary build tools if pip tries to compile (less common now). Check the official RDKit installation guide for platform-specific advice if `pip install rdkit` fails.

## 2. Runtime Errors

*   **`ValueError: numpy.dtype size changed, may indicate binary incompatibility...`**
    *   **Cause:** This is a common issue, especially in environments like Google Colab or after updating packages. It means different parts of your installed libraries were built against incompatible versions of NumPy.
    *   **Solution:**
        1.  **The primary fix:** If in Colab, **Restart the Runtime** (Runtime -> Restart runtime or Ctrl+M .). Then, **rerun all cells from the beginning**, including the installation cell.
        2.  If running locally, **recreate your virtual environment:** Deactivate (`deactivate`), remove the environment folder (`rm -rf venv`), create it again (`python -m venv venv`), activate it, and reinstall requirements (`pip install -r requirements.txt`).
        3.  Less commonly, try forcing reinstall: `pip install --upgrade --force-reinstall numpy pandas rdkit mordred`. You might still need to restart your Python kernel/environment after this.

*   **`FileNotFoundError: [Errno 2] No such file or directory: 'data/...'`**
    *   **Cause:** The script cannot find the input file specified with the `-i` argument.
    *   **Solution:**
        1.  Verify the path provided to `-i` is correct relative to where you are running the script.
        2.  Ensure the input CSV file exists at that location. Check for typos in the filename or path.

*   **`KeyError: 'Smiles'` (or your specified smiles column name)**
    *   **Cause:** The column name specified (defaulting to 'Smiles' or provided via `-s`) does not exist in the input CSV file.
    *   **Solution:**
        1.  Check the exact column names in your input CSV file (they are case-sensitive).
        2.  Use the `-s` argument to provide the correct column name if it's not `Smiles`. Example: `python calculate_descriptors.py ... -s "My SMILES Column"`

## 3. Descriptor Calculation Issues

*   **Many RDKit/Mordred Warnings (e.g., `SMILES Parse Error`, `Sanitization error`)**
    *   **Cause:** The input CSV contains invalid or problematic SMILES strings.
    *   **Solution:** This is expected for imperfect datasets. The script is designed to handle this by logging warnings, converting the problematic SMILES to `None` molecule objects, and producing `NaN` values for descriptors in the output for those rows. Review the log output and your input `data/ALLDrugswithSMILES.csv` if you need to identify and correct the specific invalid SMILES.

*   **Mordred Calculation is Very Slow**
    *   **Cause:** Mordred calculates a large number of descriptors (~1800), which can be computationally intensive, especially for large datasets.
    *   **Solution:**
        1.  Be patient. Calculation time scales with the number of molecules and the complexity of the descriptors.
        2.  Test on a smaller subset of your data first (e.g., `df.head(100).to_csv('subset.csv')`).
        3.  Ensure you are using `nproc=1` unless you are sure multiprocessing provides a benefit in your specific environment (it can sometimes add overhead or cause issues).
        4.  Run the script on a machine with sufficient RAM and CPU resources.

*   **Mordred Output Contains Many `NaN` Values**
    *   **Cause:**
        1.  The corresponding input SMILES was invalid (handled by the script).
        2.  A specific Mordred descriptor calculation failed for a valid molecule (Mordred converts these errors to `NaN`).
        3.  A descriptor might not be applicable to a particular molecule type.
        4.  If using `nproc > 1`, parallel processing errors can sometimes lead to `NaN`s.
    *   **Solution:** This is often expected. Check the log for specific Mordred errors. If entire columns are `NaN`, it might indicate a systematic issue with that descriptor calculation in your environment (see logging output in the script). Try running with `nproc=1`.

## 4. Environment Issues

*   **Inconsistent results between runs or environments:**
    *   **Cause:** Different versions of Python or the required libraries are being used.
    *   **Solution:** Always use a virtual environment and the `requirements.txt` file to ensure consistency. Check versions using `pip freeze`.

If you encounter issues not listed here, please check the script's log output for detailed error messages.
