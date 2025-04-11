import os
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

# Predefined possible bond types
all_bond_keys = ['C-C', 'C=C', 'C≡C', 'C-N', 'C=N', 'C≡N', 'C-O', 'C=O', 'C:O']

# Calculate bond length information for a molecule
def calculate_bond_lengths(mol):
    conf = mol.GetConformer()
    bond_lengths = {key: [] for key in all_bond_keys}
    
    for bond in mol.GetBonds():
        atom_idx_1 = bond.GetBeginAtomIdx()
        atom_idx_2 = bond.GetEndAtomIdx()
        bond_length = rdMolTransforms.GetBondLength(conf, atom_idx_1, atom_idx_2)

        atom1_type = mol.GetAtomWithIdx(atom_idx_1).GetSymbol()
        atom2_type = mol.GetAtomWithIdx(atom_idx_2).GetSymbol()

        if atom1_type > atom2_type:
            atom1_type, atom2_type = atom2_type, atom1_type
        
        bond_key = f'{atom1_type}-{atom2_type}'
        
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            pass
        elif bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            bond_key = f'{atom1_type}={atom2_type}'
        elif bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
            bond_key = f'{atom1_type}≡{atom2_type}'
        elif bond.GetBeginAtom().GetIsAromatic() and bond.GetEndAtom().GetIsAromatic():
            bond_key = f'{atom1_type}:{atom2_type}'
        
        if bond_key in bond_lengths:
            bond_lengths[bond_key].append(bond_length)
    
    return bond_lengths

# Process a single folder
def process(data_path_1, data_path_2, output_path, warning_log_path):
    os.makedirs(output_path, exist_ok=True)
    for folder_path, output_file_name in [(data_path_1, 'bond-length.txt'), (data_path_2, 'bond-length-optimized.txt')]:
        if not os.path.isdir(folder_path):
            print(f"Folder {folder_path} does not exist.")
            continue

        output_file_path = os.path.join(output_path, output_file_name)
        sdf_files = sorted([f for f in os.listdir(folder_path) if f.endswith('.sdf')], key=lambda x: int(x.split('.')[0]))

        with open(output_file_path, 'w') as output_file, open(warning_log_path, 'a') as warning_log:
            for sdf_file in sdf_files:
                input_sdf_file_path = os.path.join(folder_path, sdf_file)
                try:
                    suppl = Chem.SDMolSupplier(input_sdf_file_path)
                except OSError as e:
                    warning_log.write(f"Invalid input file: {input_sdf_file_path}\n")
                    continue
                
                sdf_index = sdf_file.split('.')[0]
                for mol_index, mol in enumerate(suppl):
                    if mol is not None:
                        try:
                            bond_lengths = calculate_bond_lengths(mol)
                            output_line = f"{sdf_index}\t"
                            for bond_key in all_bond_keys:
                                lengths = bond_lengths.get(bond_key, [])
                                lengths_str = ','.join(f"{length:.3f}" for length in lengths)
                                output_line += f"{bond_key}\t{lengths_str}\t"
                            output_file.write(output_line.strip() + '\n')
                        except Exception as e:
                            warning_log.write(f"Error processing {sdf_index} {mol_index}: {e}\n")
                            continue

        # Extract columns
        suffix = "-optimized" if "optimized" in output_file_name else ""
        extract_columns(output_file_path, output_path, suffix)

# Extract column information and generate new files
def extract_columns(file_path, output_dir, target_suffix=""):
    column_names = ['C-C', 'C=C', 'C≡C', 'C-N', 'C=N', 'C≡N', 'C-O', 'C=O', 'C:O']
    special_chars_map = {
        '≡': '#',
        '=': '=',
        ':': ':'
    }

    column_data = {name: [] for name in column_names}

    with open(file_path, 'r', encoding='utf-8') as file:
        for line in file:
            parts = line.strip().split('\t')
            for name in column_names:
                try:
                    index = parts.index(name) + 1
                    if index < len(parts):
                        values = parts[index].split(',')
                    else:
                        values = []
                except ValueError:
                    values = []
                column_data[name].extend(values)

    for name in column_names:
        sanitized_name = name
        for char, replacement in special_chars_map.items():
            sanitized_name = sanitized_name.replace(char, replacement)
        
        output_file = os.path.join(output_dir, f'{sanitized_name}{target_suffix}.txt')
        with open(output_file, 'w', encoding='utf-8') as file:
            non_default_values = [value for value in column_data[name] if value]
            for value in non_default_values:
                file.write(f"{value}\n")

# Main function: Traverse all folders and process
def calculate_bond_length(base_path, dataset, forcefield):
    input_folder_1 = os.path.join(base_path, dataset, 'unique_sdf')
    input_folder_2 = os.path.join(base_path, dataset, f"{forcefield}_mol")
    output_base_path = os.path.join(base_path, dataset, 'Evaluation/3D_Metrics', 'JSD/Bond')
    warning_log_path = os.path.join(base_path, dataset, 'Evaluation/3D_Metrics', 'JSD/Bond', 'warnings.log')
    if not os.path.isdir(input_folder_1) and not os.path.isdir(input_folder_2):
        print(f"No valid input folder found.")
        
    process(input_folder_1, input_folder_2, os.path.join(output_base_path), warning_log_path)

if __name__ == "__main__":
    # Set paths
    # base_path = '/data/qinyifei/model/out/'
    # output_base_path = '/data/qinyifei/model/out/3D_Metrics/JS_GEOM-Drugs/Bond/'
    base_path = '/data/qinyifei/model/teest1'
    # output_base_path = '/data/qinyifei/model/test_mol/3D_Metrics/JS_GEOM-Drugs/Bond/'
    # warning_log_path = os.path.join(output_base_path, 'warnings.log')

    # Run calculation
    calculate_bond_length(base_path, 'QM9', 'MMFF')
