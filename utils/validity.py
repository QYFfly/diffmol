import os
import csv
from rdkit import Chem

def validity_connect(input_directory, output_directory):
    # Create output directories
    os.makedirs(os.path.join(output_directory, 'QM9', 'Evaluation', 'molecular_quality'), exist_ok=True)
    os.makedirs(os.path.join(output_directory, 'QM9', 'smiles'), exist_ok=True)

    # Result output paths
    output_csv = os.path.join(output_directory, 'QM9', 'Evaluation', 'molecular_quality', 'validity_connect.csv')
    valid_smiles_output_path = os.path.join(output_directory, 'QM9', 'smiles', 'validity_largest_mol.smi')

    # Prepare CSV file
    with open(output_csv, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['Model', 'Total Molecules', 'Valid Molecules', 'Validity', 'Complete N', 'Connect'])

        # Traverse through the models in the input directory
        qm9_directory = os.path.join(input_directory)

        # Ensure it's a directory
        print(qm9_directory)
        if os.path.isdir(qm9_directory):
            # Check if the subfolder contains sdf files
            qm9_sdf_path = os.path.join(qm9_directory)

            if os.path.exists(qm9_sdf_path):
                # Process the model's SDF files and compute validity and connectivity
                validity, connect, total_molecules, valid_molecules, complete_n, valid_smiles = process_sdf_files(qm9_sdf_path, valid_smiles_output_path)

                # Write the results to CSV
                csvwriter.writerow([qm9_sdf_path, total_molecules, valid_molecules, f"{validity:.2%}", complete_n, f"{connect:.2%}"])

                # Save valid SMILES to file
                save_valid_smiles(valid_smiles_output_path, valid_smiles)

def process_sdf_files(sdf_directory, smiles_output_file):
    total_molecules = 0
    valid_molecules = 0
    complete_n = 0
    no_complete = 0
    valid_smiles = []

    # Get all SDF files
    sdf_files = [os.path.join(sdf_directory, f) for f in os.listdir(sdf_directory) if f.endswith('.sdf')]
    for sdf_file in sdf_files:
        file_index = os.path.splitext(os.path.basename(sdf_file))[0]
        supplier = Chem.SDMolSupplier(sdf_file)
        for mol in supplier:
            total_molecules += 1
            if mol is not None:
                try:
                    mol_frags = Chem.rdmolops.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
                    if len(mol_frags) > 1:
                        no_complete += 1
                    else:  # Only connected molecules are considered valid
                        largest_mol = max(mol_frags, default=mol, key=lambda m: m.GetNumAtoms())
                        Chem.SanitizeMol(largest_mol)
                        smiles = Chem.MolToSmiles(largest_mol)
                        valid_smiles.append((smiles, file_index))
                        valid_molecules += 1
                except Exception as e:
                    print(f"Error processing molecule from {sdf_file}: {e}")

    complete_n = total_molecules - no_complete
    validity = valid_molecules / total_molecules if total_molecules > 0 else 0
    connect = complete_n / total_molecules if total_molecules > 0 else 0

    return validity, connect, total_molecules, valid_molecules, complete_n, valid_smiles

def save_valid_smiles(filename, valid_smiles):
    # Ensure the save path exists
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, 'w') as f:
        for smiles, idx in valid_smiles:
            f.write(f"{smiles} {idx}\n")
