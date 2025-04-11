import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import TanimotoSimilarity
from itertools import combinations
from rdkit import rdBase

def calculate_morgan_fingerprint(molecule, radius=2, nBits=2048):
    """Calculate the Morgan fingerprint for a molecule."""
    if molecule is None:
        return None
    rdBase.DisableLog('rdApp.warning')
    fingerprint = AllChem.GetMorganFingerprintAsBitVect(molecule, radius, nBits=nBits)
    return fingerprint

def process_unique_smiles_file(unique_smiles_file):
    """Process the unique SMILES file and calculate the overall diversity based on Tanimoto similarity."""
    molecules = []
    
    with open(unique_smiles_file, 'r') as file:
        for line in file:
            smiles = line.strip().split()[0]  # Extract the SMILES string
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                molecules.append(mol)

    tanimoto_sum = 0.0
    total_molecules = len(molecules)
    total_combinations = (total_molecules * (total_molecules - 1)) / 2

    for mol1, mol2 in combinations(molecules, 2):
        fp1 = calculate_morgan_fingerprint(mol1)
        fp2 = calculate_morgan_fingerprint(mol2)
        if fp1 is None or fp2 is None:
            continue
        tanimoto_similarity = TanimotoSimilarity(fp1, fp2)
        tanimoto_sum += tanimoto_similarity
    
    overall_diversity = 1 - (tanimoto_sum / total_combinations) if total_combinations > 0 else None
    return overall_diversity

def calculate_overall_diversity(output_directory, dataset):
    """Calculate the overall diversity for a given dataset and save the result in a CSV file."""
    if dataset == 'Drug':
        smiles_dir = os.path.join(output_directory, 'Drug', 'smiles')
        output_base_folder = os.path.join(output_directory, 'Drug', 'Evaluation/molecular_quality')
    elif dataset == 'QM9':
        smiles_dir = os.path.join(output_directory, 'QM9', 'smiles')
        output_base_folder = os.path.join(output_directory, 'QM9', 'Evaluation/molecular_quality')
    else:
        print(f"Unknown dataset: {dataset}. Please specify either 'Drug' or 'QM9'.")
        return
    
    os.makedirs(output_base_folder, exist_ok=True)            
    unique_smiles_file = os.path.join(smiles_dir, 'unique_largest_mol.smi')
    output_file = os.path.join(output_base_folder, f"{dataset}_Morgan_Tdiv.csv")
    
    with open(output_file, 'w') as csvfile:
        csvfile.write('Model,Overall_Diversity\n')
        if os.path.exists(unique_smiles_file):
            model_name = os.path.basename(os.path.dirname(os.path.dirname(smiles_dir)))  # Get the model name
            overall_diversity = process_unique_smiles_file(unique_smiles_file)
            if overall_diversity is not None:
                csvfile.write(f'{model_name},{overall_diversity:.4f}\n')
            else:
                csvfile.write(f'{model_name},N/A\n')
