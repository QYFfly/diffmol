import os
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

# Predefined bond triplet types
bond_symbols_list = ['C1C-C1C-C1C', 'C12C-C12C-C12C', 'C1C-C1C-C1O', 'O1C-C1C-C1O', 'C1C-C12C-C12C', 'C1C-C2C-C1C']

def get_triple_bonds(mol):
    """Get bond triplets from the molecule"""
    valid_triple_bonds = []
    for idx_bond, bond in enumerate(mol.GetBonds()):
        idx_begin_atom = bond.GetBeginAtomIdx()
        idx_end_atom = bond.GetEndAtomIdx()
        begin_atom = mol.GetAtomWithIdx(idx_begin_atom)
        end_atom = mol.GetAtomWithIdx(idx_end_atom)
        begin_bonds = begin_atom.GetBonds()
        valid_left_bonds = [b for b in begin_bonds if b.GetIdx() != idx_bond]

        if not valid_left_bonds:
            continue

        end_bonds = end_atom.GetBonds()
        for end_bond in end_bonds:
            if end_bond.GetIdx() == idx_bond:
                continue
            for left_bond in valid_left_bonds:
                valid_triple_bonds.append([left_bond, bond, end_bond])
    return valid_triple_bonds

def find_bond_triplets(mol):
    """Calculate dihedral angles based on bond triplets"""
    bonds_list = get_triple_bonds(mol)
    angles_dict = {sym: [] for sym in bond_symbols_list}

    for bonds in bonds_list:
        sym = '-'.join([get_bond_symbol(b) for b in bonds])
        sym1 = '-'.join([get_bond_symbol(b) for b in bonds][::-1])

        if sym in angles_dict or sym1 in angles_dict:
            if sym1 in angles_dict:
                bonds = bonds[::-1]
                sym = sym1

            bond0 = bonds[0]
            atom0 = bond0.GetBeginAtomIdx()
            atom1 = bond0.GetEndAtomIdx()

            bond1 = bonds[1]
            atom1_0 = bond1.GetBeginAtomIdx()
            atom1_1 = bond1.GetEndAtomIdx()
            if atom0 == atom1_0:
                i, j, k = atom1, atom0, atom1_1
            elif atom0 == atom1_1:
                i, j, k = atom1, atom0, atom1_0
            elif atom1 == atom1_0:
                i, j, k = atom0, atom1, atom1_1
            elif atom1 == atom1_1:
                i, j, k = atom0, atom1, atom1_0
            else:
                continue

            bond2 = bonds[2]
            atom2_0 = bond2.GetBeginAtomIdx()
            atom2_1 = bond2.GetEndAtomIdx()
            if atom2_0 == k:
                l = atom2_1
            elif atom2_1 == k:
                l = atom2_0
            else:
                continue

            angle = rdMolTransforms.GetDihedralDeg(mol.GetConformer(), i, j, k, l)
            angles_dict[sym].append(angle)

    return angles_dict

def get_bond_symbol(bond):
    """Get symbolic representation of the bond"""
    a0 = bond.GetBeginAtom().GetSymbol()
    a1 = bond.GetEndAtom().GetSymbol()
    b = str(int(bond.GetBondType()))  # Single: 1, Double: 2, Triple: 3, Aromatic: 12
    return ''.join([a0, b, a1])

def process_dihedral_angles(input_folder, output_file_path, warning_log_path, suffix=""):
    """Process dihedral angle calculation for a single folder"""
    if not os.path.isdir(input_folder):
        print(f"Folder {input_folder} does not exist.")
        return

    sdf_files = sorted([f for f in os.listdir(input_folder) if f.endswith('.sdf')], 
                      key=lambda x: int(x.split('.')[0]))

    with open(output_file_path, 'w') as output_file, open(warning_log_path, 'a') as warning_log:
        for sdf_file in sdf_files:
            input_sdf_file_path = os.path.join(input_folder, sdf_file)
            try:
                suppl = Chem.SDMolSupplier(input_sdf_file_path)
            except OSError as e:
                warning_log.write(f"Invalid input file: {input_sdf_file_path}\n")
                continue
            
            sdf_index = sdf_file.split('.')[0]
            for mol_index, mol in enumerate(suppl):
                if mol is not None:
                    try:
                        angles_dict = find_bond_triplets(mol)
                        output_line = f"{sdf_index}\t"
                        for bond_sym in bond_symbols_list:
                            angles = angles_dict.get(bond_sym, [])
                            angles_str = ','.join(f"{angle:.3f}" for angle in angles)
                            output_line += f"{bond_sym}\t{angles_str}\t"
                        output_file.write(output_line.strip() + '\n')
                    except Exception as e:
                        warning_log.write(f"Error processing {sdf_index} {mol_index}: {e}\n")
                        continue

    # Extract columns
    extract_columns(output_file_path, os.path.dirname(output_file_path), suffix)

def extract_columns(file_path, output_dir, target_suffix=""):
    """Extract column info and generate new files"""
    column_data = {name: [] for name in bond_symbols_list}

    with open(file_path, 'r', encoding='utf-8') as file:
        for line in file:
            parts = line.strip().split('\t')
            for name in bond_symbols_list:
                try:
                    index = parts.index(name) + 1
                    if index < len(parts):
                        values = parts[index].split(',')
                    else:
                        values = []
                except ValueError:
                    values = []
                column_data[name].extend(values)

    for name in bond_symbols_list:
        output_file = os.path.join(output_dir, f'{name}{target_suffix}.txt')
        with open(output_file, 'w', encoding='utf-8') as file:
            non_default_values = [value for value in column_data[name] if value]
            for value in non_default_values:
                file.write(f"{value}\n")

def calculate_dihedral_angles(base_path, dataset, forcefield):
    """Main function: calculate dihedral angles"""
    # Fix forcefield path (add _mol suffix)
    forcefield_folder = f"{forcefield}_mol" if not forcefield.endswith("_mol") else forcefield
    
    input_folder_1 = os.path.join(base_path, dataset, 'unique_sdf')
    input_folder_2 = os.path.join(base_path, dataset, forcefield_folder)
    output_base_path = os.path.join(base_path, dataset, 'Evaluation/3D_Metrics/JSD/Dihedral_angles')
    warning_log_path = os.path.join(output_base_path, 'warnings.log')
    
    # Create output directory
    os.makedirs(output_base_path, exist_ok=True)
    
    # Check if input folders exist
    if not os.path.isdir(input_folder_1) and not os.path.isdir(input_folder_2):
        print(f"No valid input folders found.")
        return
    
    # Process original and optimized structures
    process_dihedral_angles(input_folder_1, 
                          os.path.join(output_base_path, 'dihedral-angles.txt'), 
                          warning_log_path)
    
    process_dihedral_angles(input_folder_2, 
                          os.path.join(output_base_path, 'dihedral-angles-optimized.txt'), 
                          warning_log_path, 
                          "-optimized")

if __name__ == "__main__":
    # Set path
    base_path = '/data/qinyifei/model/teest1'
    # Run computation
    calculate_dihedral_angles(base_path, 'QM9', 'MMFF')
