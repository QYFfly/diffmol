import os
import csv
from rdkit import Chem

class MolecularMetrics:
    def __init__(self):
        self.valid_smiles = []
        self.total_molecules = 0
        self.valid_molecules = 0
        self.complete_n = 0
        self.no_complete = 0

    def process_and_save(self, input_directory, output_directory):
        # Create output directory (including directories to save .smiles and .csv files)
        output_directory = os.path.join(output_directory, 'Drug')
        os.makedirs(output_directory, exist_ok=True)
        output_csv_dir = os.path.join(output_directory, 'Evaluation', 'molecular_quality')
        output_csv = os.path.join(output_directory, 'Evaluation', 'molecular_quality', 'validity_connect_H.csv')
        output_smiles_dir = os.path.join(output_directory, 'smiles')
        os.makedirs(output_smiles_dir, exist_ok=True)
        os.makedirs(output_csv_dir, exist_ok=True)
        
        # Prepare CSV file and write headers
        with open(output_csv, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(['Model', 'Total Molecules', 'Valid Molecules', 'Validity', 'Complete N', 'Connect'])

            # Traverse through each subfolder in the Drug folder
            drug_directory = os.path.join(input_directory)

            # Ensure it's a directory
            print(drug_directory)
            if os.path.isdir(drug_directory):
                # Check if there are sdf files in the subfolder
                drug_sdf_path = os.path.join(drug_directory)
                if os.path.exists(drug_sdf_path) and os.listdir(drug_sdf_path):  # If directory exists and has files
                    smiles_output_path = os.path.join(output_smiles_dir, f'validity_largest_mol.smi')
                    validity, connect = self.process_sdf_files(drug_sdf_path, smiles_output_path)

                    # Write results to CSV
                    csvwriter.writerow([drug_directory, self.total_molecules, self.valid_molecules, f"{validity:.2%}", self.complete_n, f"{connect:.2%}"])

    def process_sdf_files(self, sdf_directory, output_file):
        sdf_files = [os.path.join(sdf_directory, f) for f in os.listdir(sdf_directory) if f.endswith('.sdf')]
        for sdf_file in sdf_files:
            file_index = os.path.splitext(os.path.basename(sdf_file))[0]  # Extract file index
            self.process_sdf_file(sdf_file, file_index)

        # Save valid SMILES
        self.save_valid_smiles(output_file)
        self.complete_n = self.total_molecules - self.no_complete
        # Calculate validity
        validity = self.valid_molecules / self.total_molecules if self.total_molecules > 0 else 0
        connected = self.complete_n / self.total_molecules if self.total_molecules > 0 else 0
        return validity, connected

    def process_sdf_file(self, sdf_file, file_index):
        # First pass: Check molecular validity and retain hydrogen atoms
        supplier_with_hs = Chem.SDMolSupplier(sdf_file, removeHs=False)
        valid_mol_indices = []
        for i, mol in enumerate(supplier_with_hs):
            self.total_molecules += 1
            if mol is not None:
                try:
                    # Get molecular fragments
                    mol_frags = Chem.rdmolops.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
                    if len(mol_frags) > 1:
                        self.no_complete += 1
                    
                    valid_mol_indices.append(i)  # Save valid molecule indices
                except Exception as e:
                    print(f"Error processing molecule from {sdf_file}: {e}")

        # Second pass: Convert valid molecules to SMILES
        supplier = Chem.SDMolSupplier(sdf_file)
        for i in valid_mol_indices:
            mol = supplier[i]
            if mol is not None:
                try:
                    # Get molecular fragments
                    mol_frags = Chem.rdmolops.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
                    if mol_frags:
                        # Select the largest fragment
                        largest_frag = max(mol_frags, default=mol, key=lambda m: m.GetNumAtoms())
                        
                        # Sanitize the molecule and convert it to SMILES
                        Chem.SanitizeMol(largest_frag)
                        smiles = Chem.MolToSmiles(largest_frag)
                        # Save valid SMILES
                        self.valid_smiles.append((smiles, file_index))
                        self.valid_molecules += 1
                except Exception as e:
                    print(f"Error processing molecule to SMILES from {sdf_file}: {e}")

    def save_valid_smiles(self, filename):
        # Ensure the directory exists
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        with open(filename, 'w') as f:
            for smiles, idx in self.valid_smiles:
                f.write(f"{smiles} {idx}\n")
