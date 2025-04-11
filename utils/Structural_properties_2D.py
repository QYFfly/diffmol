import os
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from pathlib import Path

OUTPUT_FILE_NAME = '2D_properties.txt'

def structural_properties(file_name, smiles, output_folder):
    """
    Calculate molecular topological properties and save the results to an output file.

    :param file_name: Name of the SDF file
    :param smiles: SMILES string
    :param output_folder: Output folder path
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Unable to load molecule file: {file_name}")
            return

        # Calculate various topological properties
        num_atoms = mol.GetNumAtoms()
        num_non_hydrogen_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1)
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        num_chiral_carbons = sum(1 for idx, _ in chiral_centers if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        num_rings = rdMolDescriptors.CalcNumRings(mol)
        num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
        num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        fraction_csp3 = rdMolDescriptors.CalcFractionCSP3(mol)

        # Format the properties string
        properties_str = (f"{file_name}\tHeavy atoms:\t{num_non_hydrogen_atoms}\tChiral centers:\t{num_chiral_carbons}\t"
                          f"Number of rings:\t{num_rings}\tNumber of aromatic rings:\t{num_aromatic_rings}\t"
                          f"Number of rotatable bonds:\t{num_rotatable_bonds}\tFsp3:\t{fraction_csp3:.3f}")
        # Write the properties string to the output file
        output_path = os.path.join(output_folder, OUTPUT_FILE_NAME)
        with open(output_path, 'a') as result_file:
            result_file.write(properties_str + '\n')
    except Exception as e:
        print(f"Error processing molecule file {file_name}: {e}")

def process(input_directory, output_directory, dataset):
    """
    Process all .sdf files in the specified folder.

    :param input_directory: Root folder path
    :param folder_name: The subfolder name to process
    :param output_directory: Output folder path
    :param dataset: Dataset name ('Drug' or 'QM9')
    """
    # Dynamically select folder path
    if dataset == 'Drug':
        folder_path = os.path.join(input_directory, 'Drug', 'unique_sdf')
    elif dataset == 'QM9':
        folder_path = os.path.join(input_directory, 'QM9', 'unique_sdf')
    else:
        print(f"Unknown dataset: {dataset}. Please specify either 'Drug' or 'QM9'.")
        return

    if not os.path.isdir(folder_path):
        print(f"Folder {folder_path} does not exist, skipping.")
        return

    output_folder = os.path.join(output_directory)
    os.makedirs(output_folder, exist_ok=True)

    for file_name in os.listdir(folder_path):
        if file_name.endswith(".sdf"):
            file_path = os.path.join(folder_path, file_name)
            try:
                suppl = Chem.SDMolSupplier(file_path)
                for idx, mol in enumerate(suppl):
                    if mol is not None:
                        smiles = Chem.MolToSmiles(mol)
                        structural_properties(file_name, smiles, output_folder)
                    else:
                        print(f"Unable to parse molecule {idx} in file {file_name}")
            except Exception as e:
                print(f"Error processing file {file_name}: {e}")

def calculate_structural_properties_2D(output_directory, dataset):
    """
    Calculate the topological properties of molecules in multiple subfolders.

    :param input_directory: Root path of the data
    :param output_directory: Output path
    :param dataset: Dataset name ('Drug' or 'QM9')
    """
    overall_output_path = os.path.join(output_directory, dataset, 'Evaluation', '2D_Metrics', 'Structural_properties_2D')
    os.makedirs(overall_output_path, exist_ok=True)

    if os.path.isdir(output_directory):  # Check if folder exists
        process(output_directory, overall_output_path, dataset)
    else:
        print(f"Input folder {output_directory} does not exist, skipping.")

if __name__ == "__main__":
    input_directory = '/data/qinyifei/model/test_mol'  # Replace with actual input path
    output_directory = '/data/qinyifei/model/test_mol'  # Output path
    dataset = 'Drug'  # Can choose 'Drug' or 'QM9'

    calculate_structural_properties_2D(input_directory, output_directory, dataset)
