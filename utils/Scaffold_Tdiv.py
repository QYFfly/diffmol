import os
from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
from itertools import combinations
from rdkit.Chem import AllChem
from rdkit import rdBase
from pathlib import Path

def calculate_scaffold_similarity(smiles):
    """Extract the molecule's basic scaffold and compute its fingerprint."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetRingInfo().NumRings() == 0 or mol.GetNumAtoms() < 3:
        return None

    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None:
        return None

    scaffold_mol = Chem.MolFromSmiles(Chem.MolToSmiles(scaffold))
    rdBase.DisableLog('rdApp.warning')
    fingerprint = AllChem.GetMorganFingerprintAsBitVect(scaffold_mol, radius=2, nBits=2048)

    return fingerprint

def process_unique_smiles_file(unique_smiles_file):
    valid_smiles_list = []
    with open(unique_smiles_file, 'r') as file:
        for line in file:
            smiles = line.strip().split()[0]  # Extract SMILES string
            fp = calculate_scaffold_similarity(smiles)
            if fp is not None:
                valid_smiles_list.append(smiles)

    total_molecules = len(valid_smiles_list)
    tanimoto_sum = 0.0
    non_empty_molecules = 0

    for smiles_i, smiles_j in combinations(valid_smiles_list, 2):
        fp_i = calculate_scaffold_similarity(smiles_i)
        fp_j = calculate_scaffold_similarity(smiles_j)
        if fp_i is None or fp_j is None:
            continue
        similarity_score = DataStructs.TanimotoSimilarity(fp_i, fp_j)
        tanimoto_sum += similarity_score
        non_empty_molecules += 1

    overall_diversity = 1 - (tanimoto_sum / (total_molecules * (total_molecules - 1) / 2)) if total_molecules > 1 else 0.0
    return overall_diversity

def calculate_scaffold_diversity(output_directory, dataset):
    # Determine the path based on dataset
    if dataset == 'Drug':
        smiles_dir = Path(output_directory) / 'Drug' / 'smiles'
        output_base_folder = Path(output_directory) / 'Drug' / 'Evaluation/molecular_quality'
    elif dataset == 'QM9':
        smiles_dir = Path(output_directory) / 'QM9' / 'smiles'
        output_base_folder = Path(output_directory) / 'QM9' / 'Evaluation/molecular_quality'
    else:
        print(f"Unknown dataset: {dataset}. Please specify either 'Drug' or 'QM9'.")
        return

    os.makedirs(output_base_folder, exist_ok=True)

    unique_smiles_file = smiles_dir / 'unique_largest_mol.smi'
    output_file = output_base_folder / f"{dataset}_Scaffold_Tdiv.csv"

    with open(output_file, 'w') as csvfile:
        csvfile.write('Model,Scaffold_Diversity\n')
        if os.path.exists(unique_smiles_file):
            model_name = os.path.basename(smiles_dir.parent.parent)  # Get the model name
            scaffold_diversity = process_unique_smiles_file(unique_smiles_file)
            csvfile.write(f'{model_name},{scaffold_diversity:.4f}\n')
