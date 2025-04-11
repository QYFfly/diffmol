import os
import codecs
from rdkit.Chem import AllChem, SDMolSupplier
from rdkit import Chem
import re

def GetAngleDeg(conf, atom1, atom2, atom3):
    pt1 = conf.GetAtomPosition(atom1)
    pt2 = conf.GetAtomPosition(atom2)
    pt3 = conf.GetAtomPosition(atom3)
    v1 = [pt1.x - pt2.x, pt1.y - pt2.y, pt1.z - pt2.z]
    v2 = [pt3.x - pt2.x, pt3.y - pt2.y, pt3.z - pt2.z]
    dot_product = sum(a * b for a, b in zip(v1, v2))
    mag_v1 = sum(a**2 for a in v1)**0.5
    mag_v2 = sum(a**2 for a in v2)**0.5
    cos_theta = dot_product / (mag_v1 * mag_v2)
    angle_rad = abs(Chem.rdMolTransforms.GetAngleRad(conf, atom1, atom2, atom3))
    return angle_rad * 180 / 3.141592653589793

def get_mol_id(filename):
    match = re.search(r'([O\d+-]?\d+)\.sdf', filename)
    if match:
        return match.group(1)
    else:
        return None

def get_bond_angle_dict(mol, bond_angles_list):
    angle_dict = {}
    for bond_angles in bond_angles_list:
        substructure = Chem.MolFromSmiles(bond_angles)
        if substructure is not None:
            bond_pairs = mol.GetSubstructMatches(substructure)
            for pair in bond_pairs:
                try:
                    angle = GetAngleDeg(mol.GetConformer(), *pair)
                    assert mol.GetBondBetweenAtoms(pair[0], pair[1]) is not None
                    assert mol.GetBondBetweenAtoms(pair[2], pair[1]) is not None
                    if bond_angles not in angle_dict:
                        angle_dict[bond_angles] = []
                    angle_dict[bond_angles].append(angle)
                except Exception as e:
                    continue
    return angle_dict

def extract_columns(file_path, output_dir, target_suffix=""):
    column_names = ['CCC:', 'CSC:', 'CCO:', 'CNC:', 'OPO:', 'NCC:', 'CC=O:', 'COC:', 'CC=C:', 'OC=O:', 'NC=O:', 'CN=C:']
    
    column_data = {name: [] for name in column_names}

    with codecs.open(file_path, 'r', encoding='utf-8') as file:
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
                cleaned_values = [value.strip() for value in values if value.strip() != 'NA']
                column_data[name].extend(cleaned_values)

    os.makedirs(output_dir, exist_ok=True)
    for name in column_names:
        output_file_name = f"{name.strip(':')}{target_suffix}.txt"
        output_file_path = os.path.join(output_dir, output_file_name)
        with codecs.open(output_file_path, 'w', encoding='utf-8') as file:
            for value in column_data[name]:
                file.write(f"{value}\n")

def calculate_bond_angles(base_path, dataset, forcefield):
    bond_angles_list = ['CCC', 'CSC', 'CCO', 'CNC', 'NCC', 'CC=O', 'COC', 'CC=C', 'OC=O', 'NC=O', 'CN=C']
    
    # Adjust the forcefield path (add _mol suffix)
    forcefield_folder = f"{forcefield}_mol" if not forcefield.endswith("_mol") else forcefield
    
    input_folder_1 = os.path.join(base_path, dataset, 'unique_sdf')
    input_folder_2 = os.path.join(base_path, dataset, forcefield_folder)
    output_base_path = os.path.join(base_path, dataset, 'Evaluation/3D_Metrics/JSD/Angle')
    warning_log_path = os.path.join(output_base_path, 'warnings.log')
    
    # Create output directory
    os.makedirs(output_base_path, exist_ok=True)
    
    # Check if the input folders exist
    if not os.path.isdir(input_folder_1) and not os.path.isdir(input_folder_2):
        print(f"Valid input folder not found.")
        return
    
    def process(input_folder, output_file_path, suffix=""):
        sdf_files = [f for f in os.listdir(input_folder) if f.endswith(".sdf")]
        sorted_sdf_files = sorted(sdf_files, key=lambda x: int(get_mol_id(x)))

        with open(output_file_path, 'w', encoding='utf-8') as output_file:
            for filename in sorted_sdf_files:
                file_path = os.path.join(input_folder, filename)
                mol_id = get_mol_id(filename)
                try:
                    suppl = SDMolSupplier(file_path)
                    for mol in suppl:
                        if mol is not None:
                            try:
                                mol = Chem.AddHs(mol)
                                angle_dict = get_bond_angle_dict(mol, bond_angles_list)
                                output_file.write(f"{mol_id}\t")
                                for bond_angles in bond_angles_list:
                                    angles = angle_dict.get(bond_angles, ['NA'])
                                    output_file.write(f"{bond_angles}:\t{','.join(map(str, angles))}\t")
                                output_file.write('\n')
                            except Exception as e:
                                continue
                except OSError as e:
                    if "Invalid input file" in str(e):
                        with open(warning_log_path, 'a') as log_file:
                            log_file.write(f"{file_path}: {str(e)}\n")
                    else:
                        raise

        extract_columns(output_file_path, output_base_path, suffix)

    # Process both original and optimized structures
    process(input_folder_1, os.path.join(output_base_path, 'bond-angles.txt'))
    process(input_folder_2, os.path.join(output_base_path, 'bond-angles-optimized.txt'), "-optimized")

if __name__ == "__main__":
    # Set path
    base_path = '/data/qinyifei/model/teest1'
    # Run calculation
    calculate_bond_angles(base_path, 'QM9', 'MMFF')
