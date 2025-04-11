import os
import subprocess
import json
from multiprocessing import Pool, Manager
from pathlib import Path

def load_config(config_path):
    with open(config_path, 'r') as f:
        return json.load(f)

def delete_log_files(directory):
    """Delete log files starting with 'localhost-' or ending with '.log'."""
    log_files = [f for f in os.listdir(directory) if f.startswith("localhost-") or f.endswith(".log")]
    for log_file in log_files:
        os.remove(os.path.join(directory, log_file))

def process_sdf_file(input_sdf_file, output_sdf_file, counters, ligprep_executable):
    """Process the SDF file using the LigPrep command."""
    ligprep_command = f"{ligprep_executable} -isd {input_sdf_file} -bff 16 -g -i 0 -nt -s 1 -NJOBS 12 -osd {output_sdf_file} -NOJOBID"
    
    result = subprocess.run(ligprep_command, shell=True, 
                          stdout=subprocess.DEVNULL, 
                          stderr=subprocess.DEVNULL, 
                          text=True)
    
    counters['total'] += 1
    if result.returncode == 0 and os.path.exists(output_sdf_file):
        counters['success'] += 1

def calculate_OPLS3_preparation(output_dir: Path, dataset: str, ligprep_executable: str):
    """Prepare the OPLS3 for molecules by processing SDF files."""
    # Set up Schrödinger environment variables
    log_dir = output_dir / "logs"
    log_dir.mkdir(exist_ok=True)
    os.environ['SCHRODINGER_TMPDIR'] = str(log_dir)
    os.environ['SCHRODINGER_TEMP'] = str(log_dir)
    
    manager = Manager()
    counters = manager.dict({'total': 0, 'success': 0})

    unique_sdf_folder = output_dir / dataset / 'unique_sdf'
    if not unique_sdf_folder.exists():
        print(f"Input folder not found: {unique_sdf_folder}")
        return

    output_sdf_dir = output_dir / dataset / 'OPLS3_mol'
    output_sdf_dir.mkdir(exist_ok=True)

    sdf_files = [f for f in os.listdir(unique_sdf_folder) if f.endswith(".sdf")]
    if not sdf_files:
        print("No SDF files found in input directory")
        return

    # Process files in parallel
    with Pool(processes=os.cpu_count()) as pool:
        pool.starmap(process_sdf_file, [
            (str(unique_sdf_folder / sdf_file), 
             str(output_sdf_dir / sdf_file), 
             counters,
             ligprep_executable)
            for sdf_file in sdf_files
        ])

    # Clean up log files
    delete_log_files(str(log_dir))

    print(f"\nOPLS3 Preparation Summary for {dataset}:")
    print(f"Total files processed: {counters['total']}")
    print(f"Successfully processed files: {counters['success']}")
    print(f"Output saved to: {output_sdf_dir}")

if __name__ == "__main__":
    # Example usage
    config = load_config('/path/to/config.json')
    calculator = QM9MetricsCalculator(
        input_dir=Path(config['input_directory']),
        output_dir=Path(config['output_directory'])
    )
    
    if calculator.forcefield == "OPLS3":
        calculate_OPLS3_preparation(
            output_dir=calculator.output_dir,
            dataset=DATASET_NAME,
            ligprep_executable=calculator.ligprep_path
        )