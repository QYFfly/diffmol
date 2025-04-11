import os
import json
from multiprocessing import Pool, Manager
from rdkit import Chem
from rdkit.Chem import AllChem

def load_config(config_path):
    """Load configuration from JSON file."""
    with open(config_path, 'r') as f:
        return json.load(f)

def process_sdf_file_for_MMFF(sdf_input, sdf_output_dir, counters):
    """Process and optimize molecules in SDF files using MMFF."""
    supplier = Chem.SDMolSupplier(sdf_input)
    if not supplier:
        return

    total_molecules = 0
    embedded_molecules = 0
    optimized_molecules = 0
    output_sdf_path = os.path.join(sdf_output_dir, os.path.basename(sdf_input))
    writer = Chem.SDWriter(output_sdf_path)
    
    for mol in supplier:
        if mol is None:
            continue

        total_molecules += 1
        mol = Chem.AddHs(mol)
        try:
            embed_result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            if embed_result == -1:
                continue
        except Exception:
            continue

        embedded_molecules += 1
        try:
            optimize_result = AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94", maxIters=200)
            if optimize_result != 0:
                continue
        except Exception:
            continue

        optimized_molecules += 1
        writer.write(mol)
    
    writer.close()
    if optimized_molecules == 0 and os.path.exists(output_sdf_path):
        os.remove(output_sdf_path)
    
    with counters['lock']:
        counters['total'] += total_molecules
        counters['embedded'] += embedded_molecules
        counters['optimized'] += optimized_molecules

def calculate_MMFF_preparation(base_folder, dataset):
    """Prepare molecules in SDF files using MMFF in parallel."""
    manager = Manager()
    counters = manager.dict({'total': 0, 'embedded': 0, 'optimized': 0, 'lock': manager.Lock()})
    
    unique_sdf_folder = os.path.join(base_folder, dataset, 'unique_sdf')
    if os.path.exists(unique_sdf_folder):
        output_sdf_dir = os.path.join(base_folder, dataset, 'MMFF_mol')
        os.makedirs(output_sdf_dir, exist_ok=True)

        sdf_files = [f for f in os.listdir(unique_sdf_folder) if f.endswith(".sdf")]
        if not sdf_files:
            print("No SDF files found.")
            return

        with Pool(processes=min(40, os.cpu_count())) as pool:
            pool.starmap(
                process_sdf_file_for_MMFF,
                [(os.path.join(unique_sdf_folder, sdf_file), output_sdf_dir, counters) for sdf_file in sdf_files]
            )
    
    print(f"Total molecules: {counters['total']}")
    print(f"Embedded molecules: {counters['embedded']}")
    print(f"Optimized molecules: {counters['optimized']}")

if __name__ == "__main__":
    base_folder = '/data/qinyifei/model/teest1'
    calculate_MMFF_preparation(base_folder, 'QM9')
