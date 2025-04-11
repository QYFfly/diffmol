import os
from guacamol.distribution_matching_generator import DistributionMatchingGenerator
from guacamol.assess_distribution_learning import _assess_distribution_learning

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
# Custom generator class
class PreGeneratedMoleculeGenerator(DistributionMatchingGenerator):
    def __init__(self, smiles_file):
        """
        Initialize the generator class and load molecules from the SMILES file.
        
        :param smiles_file: Path to the SMILES file
        """
        self.smiles_file = smiles_file
        self.smiles_list = self.load_smiles()

    def load_smiles(self):
        """
        Load molecules from the SMILES file and return a list of SMILES.
        
        :return: List of SMILES molecules
        """
        with open(self.smiles_file, 'r') as f:
            smiles = [line.strip().split()[0] for line in f]
        return smiles

    def generate(self, number_samples):
        """
        Generate a given number of molecules.

        :param number_samples: Number of molecules to generate
        :return: Return the list of SMILES
        """
        return self.smiles_list

def assess_distribution_learning_with_custom_samples(generator, training_set_file, output_file, number_samples):
    """
    Evaluate molecular distribution with custom samples.

    :param generator: Molecule generator instance
    :param training_set_file: Path to the training set file
    :param output_file: Path to the output file for results
    :param number_samples: Number of samples
    """
    # Call internal evaluation function, passing the number of samples
    _assess_distribution_learning(
        model=generator,
        chembl_training_file=training_set_file,
        json_output_file=output_file,
        benchmark_version='v1',
        number_samples=number_samples
    )
    print(f"Evaluation complete, results saved to {output_file}.")

# test_file="/data/qinyifei/model/Diff-3D/data/test.smi"
# generato_test = PreGeneratedMoleculeGenerator(smiles_file=test_file)

def calculate_FCD(training_set_file, output_directory, dataset):
    """
    Traverse the directory and evaluate all molecules.

    :param base_directory: Base directory containing subdirectories for different models
    :param output_directory: Directory for output results
    :param training_set_file: Path to the training set file
    """
    input_directory = "/data/qinyifei/model/Diff-3D/data/"
    # output_directory= "/data/qinyifei/model/Diff-3D/output/"
    os.makedirs(output_directory, exist_ok=True)
    
    # if dataset == 'Drug':
    #     smiles_file = os.path.join(output_directory, 'Drug','smiles','unique_largest_mol.smi')
    #     # training_set_file = os.path.join(input_directory,'training_set_reference','Drug', 'smiles.smi')
    # elif dataset == 'QM9':
    #     smiles_file = os.path.join(output_directory, 'QM9','smiles','unique_largest_mol.smi')
    #     # training_set_file = os.path.join(input_directory,'training_set_reference','QM9', 'smiles.smi')
        
    # else:
    #     print(f"Unknown dataset: {dataset}")   
    smiles_file = os.path.join(output_directory, dataset, 'smiles', 'unique_largest_mol.smi')
    os.makedirs(output_directory, exist_ok=True)
    print(smiles_file)
    
    # Ensure it's a directory and contains the unique_largest_mol.smi file
    # test_file="/data/qinyifei/model/Diff-3D/data/test.smi"
    # generato_test = PreGeneratedMoleculeGenerator(smiles_file=test_file)
    # if os.path.isdir(smiles_file):
        # Create the custom generator instance
    generator = PreGeneratedMoleculeGenerator(smiles_file=smiles_file)

    # Calculate the number of molecules in the generated molecule file
    num_molecules = len(generator.smiles_list)
    print(f"Processing, number of generated molecules: {num_molecules}")

    # Dynamically generate the output file path
    output_file = os.path.join(output_directory, dataset, 'Evaluation', 'Distribution_similarity_score', f'{dataset}FCD.json')
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Evaluate the generated molecules, using the number of molecules in the file
    assess_distribution_learning_with_custom_samples(generator, training_set_file, output_file, num_molecules)

# input_directory= "/data/qinyifei/model/Diff-3D/data/"
# output_directory= "/data/qinyifei/model/Diff-3D/output/"
# calculate_FCD(input_directory, output_directory, 'QM9')
