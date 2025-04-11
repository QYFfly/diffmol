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

    # Calculate KDE and add smoothing factor
    kde_P = gaussian_kde(data1)
    kde_Q = gaussian_kde(data2)
    
    # Create evaluation range (extend 5% beyond the boundaries)
    min_val = min(data1.min(), data2.min())
    max_val = max(data1.max(), data2.max())
    range_extension = (max_val - min_val) * 0.05
    x_eval = np.linspace(min_val - range_extension, 
                        max_val + range_extension, 
                        num=num_points)
    
    # Calculate probability distributions
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
    data = data[~np.isnan(data)]  # Remove NaN values
    
    if len(data) == 0:
        raise ValueError("The file contains no valid data")
    
    return data

def calculate_JSD_dihedral_angles(base_path, dataset):
    """Calculate JS Divergence for Dihedral Angles"""
    # Define file pairs to compare
    file_pairs = [
        ('C1C-C1C-C1C-optimized.txt', 'C1C-C1C-C1C.txt'),
        ('C12C-C12C-C12C-optimized.txt', 'C12C-C12C-C12C.txt'),
        ('C1C-C1C-C1O-optimized.txt', 'C1C-C1C-C1O.txt'),
        ('O1C-C1C-C1O-optimized.txt', 'O1C-C1C-C1O.txt'),
        ('C1C-C12C-C12C-optimized.txt', 'C1C-C12C-C12C.txt'),
        ('C1C-C2C-C1C-optimized.txt', 'C1C-C2C-C1C.txt')
    ]

    # Set paths
    input_dir = os.path.join(base_path, dataset, 'Evaluation/3D_Metrics/JSD/Dihedral_angles')
    output_dir = os.path.join(base_path, dataset, 'Evaluation/3D_Metrics/JSD/results')
    os.makedirs(output_dir, exist_ok=True)
    
    # Result file paths
    result_file = os.path.join(output_dir, 'dihedral_angle_jsd_results.txt')
    error_log = os.path.join(output_dir, 'dihedral_error_log.txt')
    
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
            dihedral_type = pair[1].replace('.txt', '')
            results.append(f"{dihedral_type}\t{jsd:.4f}\n")
            
        except Exception as e:
            errors.append(f"{pair[1]}: {str(e)}\n")
            continue

    # Write results to file
    with open(result_file, 'w') as f:
        f.write("DihedralType\tJSD_Value\n")
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
    calculate_JSD_dihedral_angles(base_path, dataset)
