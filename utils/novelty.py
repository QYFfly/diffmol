import os
import csv
import json
class MolecularMetrics:
    def __init__(self, train_smiles_file, unique_smiles_file):
        self.train_smiles_file = train_smiles_file
        self.unique_smiles_file = unique_smiles_file
        self.train_smiles = self.load_smiles(self.train_smiles_file)

    def load_smiles(self, filepath, has_index=False):
        """Load the SMILES from the file."""
        smiles_list = []
        with open(filepath, 'r') as file:
            for line in file:
                if has_index:
                    smiles = line.strip().split()[0]  # Extract the first part as SMILES
                else:
                    smiles = line.strip()
                smiles_list.append(smiles)
        return smiles_list

    def compute_novelty(self):
        """Compute the novelty of the dataset."""
        num_novel = 0
        novel = []
        unique_smiles = self.load_smiles(self.unique_smiles_file, has_index=True)

        if self.train_smiles is None:
            print("Dataset smiles is None, novelty computation skipped")
            return 1, 1

        for smiles in unique_smiles:
            if smiles not in self.train_smiles:
                novel.append(smiles)
                num_novel += 1
        
        novelty_ratio = num_novel / len(unique_smiles) if len(unique_smiles) > 0 else 0
        # print(novelty_ratio)
        return num_novel, len(unique_smiles), novelty_ratio

def computely_novelty(output_directory, dataset, train_smiles_file):
    """
    Calculate the novelty ratio and save it as a CSV file
    :param input_directory: Input file directory
    :param output_directory: Output file directory
    :param dataset: Dataset type (QM9 or Drug)
    """
    if dataset == 'QM9':
        # train_smiles_file = '/data/qinyifei/model/Diff-3D/data/data_reference/QM9/QM9_smiles.txt'
        unique_smiles_file = os.path.join(output_directory, 'QM9', 'smiles', 'unique_largest_mol.smi')
        output_csv = os.path.join(output_directory, 'QM9', 'Evaluation','molecular_quality', 'QM9_novelty_rate.csv')
    elif dataset == 'Drug':
        # train_smiles_file = '/data/qinyifei/model/Diff-3D/data/data_reference/Drug/Drug_smiles.smi'
        unique_smiles_file = os.path.join(output_directory, 'Drug', 'smiles', 'unique_largest_mol.smi')
        output_csv = os.path.join(output_directory, 'Drug','Evaluation','molecular_quality', 'Drug_novelty_rate.csv')
        
    else:
        print(f"Unknown dataset: {dataset}. Please specify either 'QM9' or 'Drug'.")
        return

    # Prepare CSV file for output
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    with open(output_csv, 'w', newline='') as csvfile:
        # print(output_csv)
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['Model', 'Num Novel', 'Unique Smiles Count', 'Novelty Ratio'])

        # Traverse all model folders
    
        # print(unique_smiles_file)
        if os.path.exists(unique_smiles_file):
            # print(unique_smiles_file)
            metrics = MolecularMetrics(train_smiles_file, unique_smiles_file)
            num_novel, unique_count, novelty_ratio = metrics.compute_novelty()

            # Write results to CSV
            csvwriter.writerow([num_novel, unique_count, f"{novelty_ratio:.2%}"])

def main():
    # Call compute_novelty to process the data
    # Load `dataset`, `input_directory`, and `output_directory` from the config file
    config_path = '/data/qinyifei/model/Diff-3D/config.json'
    with open(config_path, 'r') as f:
        config = json.load(f)

    input_directory = config.get('output_directory', '/data/qinyifei/model/Diff-3D/output/')
    output_directory = config.get('output_directory', '/data/qinyifei/model/Diff-3D/output/')
    dataset = config.get('dataset', 'QM9')  # Default to QM9

    # Calculate novelty
    computely_novelty(input_directory, output_directory, dataset)

if __name__ == "__main__":
    main()
