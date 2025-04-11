import os
import glob
from rdkit import Chem
from rdkit.Chem import QED
import pandas as pd

def process_sdf_file(sdf_file):
    """
    Process a single SDF file and compute the QED value for each molecule.
    """
    suppl = Chem.SDMolSupplier(sdf_file)
    results = []

    # Get the file prefix as the ID, e.g., for 6093.sdf, the ID is 6093
    file_prefix = os.path.splitext(os.path.basename(sdf_file))[0]

    for mol in suppl:
        if mol is not None:
            qed_value = QED.qed(mol)
            results.append((file_prefix, qed_value))
        else:
            print(f"Warning: Molecule in file {sdf_file} could not be parsed.")

    return results

def calculate_QED_for_subfolders(output_directory, dataset):
    """
    Calculate the QED values for the given dataset and save them to a CSV file.
    """

    # Choose unique_sdf folder path based on dataset
    if dataset == 'Drug':
        unique_sdf_folder = os.path.join(output_directory, 'Drug', 'unique_sdf')
        output_base_folder = os.path.join(output_directory, 'Drug', 'Evaluation/molecular_quality')
    elif dataset == 'QM9':
        unique_sdf_folder = os.path.join(output_directory, 'QM9', 'unique_sdf')
        output_base_folder = os.path.join(output_directory, 'QM9', 'Evaluation/molecular_quality')
    else:
        print(f"Unknown dataset: {dataset}. Please specify either 'Drug' or 'QM9'.")
        return
    os.makedirs(output_base_folder, exist_ok=True)
    
    # Check if unique_sdf folder exists
    if not os.path.exists(unique_sdf_folder):
        print(f"Error: The directory {unique_sdf_folder} does not exist.")
        return

    # Find all SDF files
    sdf_files = glob.glob(os.path.join(unique_sdf_folder, '*.sdf'))

    if not sdf_files:
        print(f"No SDF files found in {unique_sdf_folder}.")
        return

    # Create output file path
    output_file = os.path.join(output_base_folder, f"{dataset}_QED_values.csv")

    results = []
    for sdf_file in sdf_files:
        results.extend(process_sdf_file(sdf_file))

    # Save results to a CSV file
    df = pd.DataFrame(results, columns=['ID', 'QED'])
    df.to_csv(output_file, index=False)
    print(f"Saved QED values to {output_file}")
