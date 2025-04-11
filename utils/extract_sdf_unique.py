import os
import shutil
import json

def process_unique_sdf(input_directory, output_directory, dataset):
    """
    Process the unique_largest_mol.smi file in each model and copy the corresponding .sdf file to the unique_sdf folder
    :param input_directory: Input directory
    :param output_directory: Output directory
    :param dataset: Dataset type (QM9 or Drug)
    """
    # Set paths based on dataset
    if dataset == 'QM9':
        model_base_folder = os.path.join(input_directory)
        smiles_base_path = os.path.join(output_directory, 'QM9')
        output_base_folder = os.path.join(output_directory, 'QM9', 'unique_sdf')
    elif dataset == 'Drug':
        model_base_folder = os.path.join(input_directory)
        smiles_base_path = os.path.join(output_directory, 'Drug')
        output_base_folder = os.path.join(output_directory, 'Drug', 'unique_sdf')
    else:
        print(f"Unknown dataset: {dataset}. Please specify either 'QM9' or 'Drug'.")
        return

    # Ensure the output folder exists
    if not os.path.exists(output_base_folder):
        os.makedirs(output_base_folder)

    # Iterate over each model folder
    
    model_folder_path = os.path.join(model_base_folder)
    if os.path.isdir(model_folder_path):
        # Locate the unique_largest_mol.smi file
        smiles_file_path = os.path.join(smiles_base_path, 'smiles', 'unique_largest_mol.smi')

        if os.path.exists(smiles_file_path):
            with open(smiles_file_path, 'r') as smiles_file:
                lines = smiles_file.readlines()
                
                # Extract the numeric ID from each line
                numbers = [line.strip().split()[-1] for line in lines]  # Extract the last part as numeric ID
                
                # Iterate over each numeric ID, find the corresponding .sdf file and copy it to the unique_sdf folder
                for number in numbers:
                    sdf_file_path = os.path.join(model_folder_path, f'{number}.sdf')
                    if os.path.exists(sdf_file_path):
                        # Copy to output_base_folder
                        shutil.copy(sdf_file_path, output_base_folder)
                        # print(f"Copy successful: {sdf_file_path} to {output_base_folder}")
                    else:
                        print(f"Could not find .sdf file: {sdf_file_path}")
        else:
            print(f"Could not find unique_largest_mol.smi file at: {smiles_file_path}")
    else:
        print(f"Could not find model folder: {model_folder_path}")

def main():
    # Read the configuration file to get input and output paths and dataset type
    config_path = '/data/qinyifei/model/Diff-3D/config.json'
    with open(config_path, 'r') as f:
        config = json.load(f)

    input_directory = config.get('input_directory', '/data/qinyifei/model/Diff-3D/data/')
    output_directory = config.get('output_directory', '/data/qinyifei/model/Diff-3D/output/')
    dataset = config.get('dataset', 'QM9')  # Default is QM9

    # Process the unique_largest_mol.smi and copy the corresponding .sdf files
    process_unique_sdf(input_directory, output_directory, dataset)

if __name__ == "__main__":
    main()
