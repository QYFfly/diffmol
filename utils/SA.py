import os
import glob
import gzip
import math
import pickle
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import pandas as pd

_fscores = None

def readFragmentScores(forcefield_path):
    """
    Load fragment scores from a specified .pkl.gz file.
    """
    global _fscores
    forcefield_path = str(forcefield_path)  # Ensure compatibility with Path or str
    data = pickle.load(gzip.open(forcefield_path, 'rb'))
    outDict = {}
    for i in data:
        for j in range(1, len(i)):
            outDict[i[j]] = float(i[0])
    _fscores = outDict

def numBridgeheadsAndSpiro(mol, ri=None):
    nSpiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
    nBridgehead = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
    return nBridgehead, nSpiro

def calculateSA(mol, forcefield_path):
    if _fscores is None:
        readFragmentScores(forcefield_path)

    if mol is None:
        return None

    fp = rdMolDescriptors.GetMorganFingerprint(mol, 2)
    fps = fp.GetNonzeroElements()
    score1 = 0.
    nf = 0

    for bitId, v in fps.items():
        nf += v
        sfp = bitId
        score1 += _fscores.get(sfp, -4) * v
    score1 /= nf

    nAtoms = mol.GetNumAtoms()
    nChiralCenters = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    ri = mol.GetRingInfo()
    nBridgeheads, nSpiro = numBridgeheadsAndSpiro(mol, ri)
    nMacrocycles = 0
    for x in ri.AtomRings():
        if len(x) > 8:
            nMacrocycles += 1

    sizePenalty = nAtoms**1.005 - nAtoms
    stereoPenalty = math.log10(nChiralCenters + 1)
    spiroPenalty = math.log10(nSpiro + 1)
    bridgePenalty = math.log10(nBridgeheads + 1)
    macrocyclePenalty = math.log10(2) if nMacrocycles > 0 else 0.

    score2 = 0. - sizePenalty - stereoPenalty - spiroPenalty - bridgePenalty - macrocyclePenalty

    score3 = 0.
    if nAtoms > len(fps):
        score3 = math.log(float(nAtoms) / len(fps)) * .5

    sascore = score1 + score2 + score3

    min = -4.0
    max = 2.5
    sascore = 11. - (sascore - min + 1) / (max - min) * 9.

    if sascore > 8.:
        sascore = 8. + math.log(sascore + 1. - 9.)
    if sascore > 10.:
        sascore = 10.0
    elif sascore < 1.:
        sascore = 1.0

    return sascore

def process_sdf_file(sdf_file, forcefield_path):
    suppl = Chem.SDMolSupplier(sdf_file)
    results = []

    for mol in suppl:
        if mol is None:
            continue
        
        sa_score = calculateSA(mol, forcefield_path)
        if sa_score is not None:
            molecule_id = os.path.basename(sdf_file).replace('.sdf', '')
            results.append((molecule_id, sa_score))

    return results

def calculate_SA_for_subfolders(output_directory, dataset, forcefield_path):
    """
    Calculate the SA values for the given dataset and save them to a CSV file.
    """
    if dataset == 'Drug':
        unique_sdf_folder = os.path.join(output_directory, 'Drug', 'unique_sdf')
        output_base_folder = os.path.join(output_directory, 'Drug', 'Evaluation/molecular_quality')
    elif dataset == 'QM9':
        unique_sdf_folder = os.path.join(output_directory, 'QM9', 'unique_sdf')
        output_base_folder = os.path.join(output_directory, 'QM9', 'Evaluation/molecular_quality')
    else:
        print(f"Unknown dataset: {dataset}. Please specify either 'Drug' or 'QM9'.")
        return

    os.makedirs(output_base_folder, exist_ok=True)

    if not os.path.exists(unique_sdf_folder):
        print(f"Error: The directory {unique_sdf_folder} does not exist.")
        return

    sdf_files = glob.glob(os.path.join(unique_sdf_folder, '*.sdf'))

    if not sdf_files:
        print(f"No SDF files found in {unique_sdf_folder}.")
        return

    output_file = os.path.join(output_base_folder, f"{dataset}_SA_values.csv")

    results = []
    for sdf_file in sdf_files:
        results.extend(process_sdf_file(sdf_file, forcefield_path))

    df = pd.DataFrame(results, columns=['ID', 'SA_Score'])
    df.to_csv(output_file, index=False)
    print(f"Saved SA values to {output_file}")
