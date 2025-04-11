import os
import csv
import json

def remove_duplicates(smiles_file, output_file):
    unique_smiles = set()
    total_molecules = 0

    with open(smiles_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            total_molecules += 1
            smiles = line.split()[0]  # 仅提取 SMILES 字符串部分
            if smiles not in unique_smiles:
                unique_smiles.add(smiles)
                outfile.write(line)

    unique_count = len(unique_smiles)
    unique_rate = unique_count / total_molecules if total_molecules > 0 else 0

    return total_molecules, unique_count, unique_rate

def process_all_models(output_directory, dataset, output_csv):
    # Prepare CSV file for output
    if not os.path.exists(os.path.dirname(output_csv)):
        os.makedirs(os.path.dirname(output_csv))

    with open(output_csv, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['Model', 'Total Molecules', 'Unique Molecules', 'Unique Rate'])

        # Based on dataset, adjust the directory structure
        dataset_directory = os.path.join(output_directory, dataset)
        smiles_directory = os.path.join(dataset_directory, 'smiles')
        
        
           
            
        input_smiles = os.path.join(smiles_directory, 'validity_largest_mol.smi')
        output_smiles = os.path.join(smiles_directory, 'unique_largest_mol.smi')
        print(input_smiles,output_smiles)
        if os.path.exists(input_smiles):
            total_molecules, unique_count, unique_rate = remove_duplicates(input_smiles, output_smiles)

            # Write results to CSV
            csvwriter.writerow([total_molecules, unique_count, f"{unique_rate:.2%}"])

def load_config(config_path):
    with open(config_path, 'r') as f:
        config = json.load(f)
    return config

def calculate_uniqueness(output_directory, dataset):
    # Output CSV path
    output_csv = os.path.join(output_directory, dataset, 'Evaluation','molecular_quality', 'unique_rate.csv')

    # Process all models and save unique rate statistics
    process_all_models(output_directory, dataset, output_csv)
