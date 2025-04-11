import numpy as np
from scipy.spatial.distance import jensenshannon
import os
from scipy.stats import gaussian_kde

def calculate_js_divergence(data1, data2, num_points=100, epsilon=1e-6):
    """Calculate JS Divergence between two datasets"""
    if len(data1) < 2 or len(data2) < 2:
        raise ValueError("Not enough data points, at least 2 data points are required for KDE estimation")
    
    if np.all(data1 == data1[0]) or np.all(data2 == data2[0]):
        raise ValueError("Data is too homogeneous (all values are the same)")

    kde_P = gaussian_kde(data1)
    kde_Q = gaussian_kde(data2)
    x_eval = np.linspace(min(data1.min(), data2.min()), 
                        max(data1.max(), data2.max()), 
                        num=num_points)
    P = kde_P(x_eval) + epsilon
    Q = kde_Q(x_eval) + epsilon
    P /= np.sum(P)
    Q /= np.sum(Q)
    return jensenshannon(P, Q)

def read_data(file_path):
    """Read data from file and handle NaN values"""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
    
    data = np.genfromtxt(file_path)
    data = data[~np.isnan(data)]
    
    if len(data) == 0:
        raise ValueError("The file contains no valid data")
    
    return data

def calculate_JSD_angles(base_path, dataset):
    """Calculate JS divergence for bond angles"""
    # Define file pairs to compare
    file_pairs = [
        ('CCC-optimized.txt', 'CCC.txt'),
        ('CSC-optimized.txt', 'CSC.txt'),
        ('CCO-optimized.txt', 'CCO.txt'),
        ('CNC-optimized.txt', 'CNC.txt'), 
        ('NCC-optimized.txt', 'NCC.txt'),
        ('CC=O-optimized.txt', 'CC=O.txt'),
        ('COC-optimized.txt', 'COC.txt'),
        ('CC=C-optimized.txt', 'CC=C.txt'),
        ('OC=O-optimized.txt', 'OC=O.txt'),
        ('NC=O-optimized.txt', 'NC=O.txt'),
        ('CN=C-optimized.txt', 'CN=C.txt')
    ]

    # Set paths
    input_dir = os.path.join(base_path, dataset, 'Evaluation/3D_Metrics/JSD/Angle')
    output_dir = os.path.join(base_path, dataset, 'Evaluation/3D_Metrics/JSD/results')
    os.makedirs(output_dir, exist_ok=True)
    
    # Result file paths
    result_file = os.path.join(output_dir, 'bond_angles_jsd_results.txt')
    error_log = os.path.join(output_dir, 'angle_error_log.txt')
   
    results = []
    errors = []

    for pair in file_pairs:
        try:
            # Read data
            data_opt = read_data(os.path.join(input_dir, pair[0]))
            data_raw = read_data(os.path.join(input_dir, pair[1]))
            
            # Calculate JS divergence
            jsd = calculate_js_divergence(data_opt, data_raw)
            
            # Record the result
            angle_type = pair[1].replace('.txt', '')
            results.append(f"{angle_type}\t{jsd:.4f}\n")
            
        except Exception as e:
            errors.append(f"{pair[1]}: {str(e)}\n")
            continue

    # Write results to file
    with open(result_file, 'w') as f:
        f.write("AngleType\tJSD_Value\n")
        f.writelines(results)
    
    # Log errors
    if errors:
        with open(error_log, 'w') as f:
            f.writelines(errors)
        print(f"Calculation completed with {len(errors)} errors, see {error_log} for details")
    else:
        print(f"Calculation completed, results saved to {result_file}")

if __name__ == "__main__":
    # Configuration parameters
    base_path = '/data/qinyifei/model/teest1/'
    dataset = 'QM9'  # Specify the dataset name
    
    # Run the calculation
    calculate_JSD_angles(base_path, dataset)
