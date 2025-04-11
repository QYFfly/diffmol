import os

# Define the properties to extract and their corresponding filenames
properties = {
    'Heavy_atoms': 'Heavy_atoms.txt',
    'Chiral_centers': 'Chiral_centers.txt',
    'Rings': 'Rings.txt',
    'Aromatic_rings': 'Aromatic_rings.txt',
    'Rotatable_bonds': 'Rotatable_bonds.txt',
    'Fsp3': 'Fsp3.txt'
}

def extract_structural_properties(input_directory, dataset):
    """
    Extract molecular structural properties from the specified dataset (Drug or QM9) and save them to files.

    :param dataset: Dataset name, 'Drug' or 'QM9'
    :param input_directory: Root path of the data directory
    """
    # Set folder path based on dataset
    if dataset == 'Drug':
        folder_path = os.path.join(input_directory, 'Drug', 'Evaluation', '2D_Metrics', 'Structural_properties_2D')
    elif dataset == 'QM9':
        folder_path = os.path.join(input_directory, 'QM9', 'Evaluation', '2D_Metrics', 'Structural_properties_2D')
    else:
        print(f"Unknown dataset: {dataset}")
        folder_path = None

    # Ensure the folder path is valid
    if folder_path and os.path.isdir(folder_path):
        # Traverse the folder and process data
        for subdir, dirs, files in os.walk(folder_path):
            for file in files:
                if file == '2D_properties.txt':
                    file_path = os.path.join(subdir, file)
                    with open(file_path, 'r') as f:
                        lines = f.readlines()

                    # Initialize dictionary to store property data
                    data = {key: [] for key in properties.keys()}

                    # Process data line by line
                    for line in lines:
                        parts = line.strip().split('\t')
                        try:
                            data['Heavy_atoms'].append(parts[2])
                            data['Chiral_centers'].append(parts[4])
                            data['Rings'].append(parts[6])
                            data['Aromatic_rings'].append(parts[8])
                            data['Rotatable_bonds'].append(parts[10])
                            data['Fsp3'].append(parts[12])
                        except (IndexError, ValueError):
                            print(f"Error processing line: {line}")
                            continue

                    # Write data to respective files
                    for prop, filename in properties.items():
                        output_file = os.path.join(subdir, filename)
                        with open(output_file, 'w') as out_f:
                            out_f.write('\n'.join(data[prop]))

        print("Data extraction completed and saved to respective files.")
    else:
        print(f"Input path {folder_path} does not exist, skipping processing.")
