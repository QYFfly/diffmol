import os
import numpy as np
from rdkit import Chem
import csv

# Parse SDF file and extract atom coordinates and types
def parse_sdf(file_path, COUNT):
    supplier = Chem.SDMolSupplier(file_path, removeHs=False)  # Ensure hydrogen atoms are not removed
    positions_list = []
    atom_types_list = []
    for mol in supplier:
        if mol is None:
            COUNT = COUNT + 1
            continue
        positions = mol.GetConformer().GetPositions()
        atom_types = [mol.GetAtomWithIdx(i).GetSymbol() for i in range(mol.GetNumAtoms())]
        positions_list.append(positions)
        atom_types_list.append(atom_types)
    return positions_list, atom_types_list, COUNT

# Define bond length data and functions
bonds1 = {'H': {'H': 74, 'C': 109, 'N': 101, 'O': 96, 'F': 92, 'B': 119, 'Si': 148, 'P': 144, 'As': 152, 'S': 134, 'Cl': 127, 'Br': 141, 'I': 161},
          'C': {'H': 109, 'C': 154, 'N': 147, 'O': 143, 'F': 135, 'Si': 185, 'P': 184, 'S': 182, 'Cl': 177, 'Br': 194, 'I': 214},
          'N': {'H': 101, 'C': 147, 'N': 145, 'O': 140, 'F': 136, 'Cl': 175, 'Br': 214, 'S': 168, 'I': 222, 'P': 177},
          'O': {'H': 96, 'C': 143, 'N': 140, 'O': 148, 'F': 142, 'Br': 172, 'S': 151, 'P': 163, 'Si': 163, 'Cl': 164, 'I': 194},
          'F': {'H': 92, 'C': 135, 'N': 136, 'O': 142, 'F': 142, 'S': 158, 'Si': 160, 'Cl': 166, 'Br': 178, 'P': 156, 'I': 187},
          'B': {'H':  119, 'Cl': 175},
          'Si': {'Si': 233, 'H': 148, 'C': 185, 'O': 163, 'S': 200, 'F': 160, 'Cl': 202, 'Br': 215, 'I': 243 },
          'Cl': {'Cl': 199, 'H': 127, 'C': 177, 'N': 175, 'O': 164, 'P': 203, 'S': 207, 'B': 175, 'Si': 202, 'F': 166, 'Br': 214},
          'S': {'H': 134, 'C': 182, 'N': 168, 'O': 151, 'S': 204, 'F': 158, 'Cl': 207, 'Br': 225, 'Si': 200, 'P': 210, 'I': 234},
          'Br': {'Br': 228, 'H': 141, 'C': 194, 'O': 172, 'N': 214, 'Si': 215, 'S': 225, 'F': 178, 'Cl': 214, 'P': 222},
          'P': {'P': 221, 'H': 144, 'C': 184, 'O': 163, 'Cl': 203, 'S': 210, 'F': 156, 'N': 177, 'Br': 222},
          'I': {'H': 161, 'C': 214, 'Si': 243, 'N': 222, 'O': 194, 'S': 234, 'F': 187, 'I': 266},
          'As': {'H': 152}
          }

bonds2 = {'C': {'C': 134, 'N': 129, 'O': 120, 'S': 160},
          'N': {'C': 129, 'N': 125, 'O': 121},
          'O': {'C': 120, 'N': 121, 'O': 121, 'P': 150},
          'P': {'O': 150, 'S': 186},
          'S': {'P': 186}}

bonds3 = {'C': {'C': 120, 'N': 116},
          'N': {'C': 116, 'N': 110},
        }

allowed_bonds = {'H': 1, 'C': 4, 'N': 3, 'O': 2, 'F': 1, 'B': 3, 'Al': 3,
                 'Si': 4, 'P': [3, 5],
                 'S': 4, 'Cl': 1, 'As': 3, 'Br': 1, 'I': 1, 'Hg': [1, 2],
                 'Bi': [3, 5]}
margin1, margin2, margin3 = 10, 5, 3

def get_bond_order(atom1, atom2, distance, check_exists=False):
    distance = 100 * distance  # Change the metric to pm
    if check_exists:
        if atom1 not in bonds1:
            return 0
        if atom2 not in bonds1[atom1]:
            return 0

    if distance < bonds1[atom1][atom2] + margin1:
        if atom1 in bonds2 and atom2 in bonds2[atom1]:
            if distance < bonds2[atom1][atom2] + margin2:
                if atom1 in bonds3 and atom2 in bonds3[atom1]:
                    if distance < bonds3[atom1][atom2] + margin3:
                        return 3  # Triple bond
                return 2  # Double bond
        return 1  # Single bond
    return 0  # No bond

def check_stability(positions, atom_types, debug=False):
    nr_atoms = len(atom_types)
    atom_type = atom_types
    assert len(positions.shape) == 2
    assert positions.shape[1] == 3
    x = positions[:, 0]
    y = positions[:, 1]
    z = positions[:, 2]
    nr_bonds = np.zeros(len(x), dtype='int')
    for i in range(len(x)):
        for j in range(i + 1, len(x)):
            p1 = np.array([x[i], y[i], z[i]])
            p2 = np.array([x[j], y[j], z[j]])
            dist = np.sqrt(np.sum((p1 - p2) ** 2))
            
            pair = sorted([atom_type[i], atom_type[j]])  # Ensure the order is consistent
            order = get_bond_order(pair[0], pair[1], dist, check_exists=True)
            nr_bonds[i] += order
            nr_bonds[j] += order

    nr_stable_bonds = 0
    for atom_type_i, nr_bonds_i in zip(atom_type, nr_bonds):
        possible_bonds = allowed_bonds[atom_type_i]
        if type(possible_bonds) == int:
            is_stable = possible_bonds == nr_bonds_i
        else:
            is_stable = nr_bonds_i in possible_bonds
        if not is_stable and debug:
            print(f"Invalid bonds for atom {atom_type_i} with {nr_bonds_i} bonds")
        nr_stable_bonds += int(is_stable)

    molecule_stable = nr_stable_bonds == len(atom_types)
    return molecule_stable, nr_stable_bonds, len(atom_types)

def calculate_stable( output_directory, dataset, debug=False):
    # Set folder path according to dataset
    if dataset == 'Drug':
        unique_sdf_folder = os.path.join(output_directory, 'Drug', 'unique_sdf')
        output_folder = os.path.join(output_directory, 'Drug', 'Evaluation/3D_Metrics')
    elif dataset == 'QM9':
        unique_sdf_folder = os.path.join(output_directory, 'QM9', 'unique_sdf')
        output_folder= os.path.join(output_directory, 'QM9', 'Evaluation/3D_Metrics')
       
    else:
        print(f"Unknown dataset: {dataset}. Please specify either 'Drug' or 'QM9'.")
        return
    os.makedirs(output_folder, exist_ok=True)
    output_file=os.path.join(output_folder, f"{dataset}_stable.csv")
    # Open output CSV file and prepare to write results
    with open(output_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['stable_molecules', 'total_molecules', 'stability_ratio', 'stability_Atom_ratio'])
        
        sum_stable = 0  # Count stable molecules
        sum_total = 0   # Count total molecules
        sum_Atom_stable = 0  # Count stable atoms
        sum_Atom_total = 0    # Count total atoms
        COUNT = 0  # Count skipped invalid files
        
        # Iterate through each SDF file in the model folder
        for sdf_file in os.listdir(unique_sdf_folder):
            if sdf_file.endswith('.sdf'):  # Process only SDF files
                file_path = os.path.join(unique_sdf_folder, sdf_file)
                positions_list, atom_types_list, COUNT = parse_sdf(file_path, COUNT)  # Parse SDF file
                
                # Iterate through each molecule
                for positions, atom_types in zip(positions_list, atom_types_list):
                    molecule_stable, nr_stable_bonds, nr_atoms = check_stability(positions, atom_types, debug=debug)  # Check stability
                    sum_total += 1  # Increase total molecule count by 1
                    sum_Atom_stable += nr_stable_bonds  # Add stable atom count
                    sum_Atom_total += nr_atoms  # Add total atom count
                    
                    # If molecule is stable, increase stable molecule count by 1
                    if molecule_stable:
                        sum_stable += 1
        
        # Calculate stability ratio
        if sum_total > 0:
            stability_ratio = sum_stable / sum_total
        else:
            stability_ratio = 0

        # Calculate stable atom ratio
        if sum_Atom_total > 0:
            stability_Atom_ratio = sum_Atom_stable / sum_Atom_total
        else:
            stability_Atom_ratio = 0

        # Write results to CSV file
        csv_writer.writerow([ sum_stable, sum_total, stability_ratio, stability_Atom_ratio])

        # Print processing information
        print(f"Processed : Stable molecules {sum_stable}/{sum_total}, Stability ratio: {stability_ratio}, "
                f"Stable Atom {sum_Atom_stable}/{sum_Atom_total}, stability_Atom_ratio: {stability_Atom_ratio}")


if __name__ == "__main__":
    # Input root directory and output file path
    root_dir = "/data/qinyifei/model/test_mol"
    output_file = "/data/qinyifei/model/test_mol/2D_Metrics/GEOM-Drugs_stable.csv"
    
    # Call the function to process
    calculate_stable(root_dir, output_file, debug=False)
