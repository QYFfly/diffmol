# diffmol
DiffMol is a diffusion-model-based 3D molecular generation and evaluation framework that enables quantitative assessment of key metrics including geometric rationality, chemical properties, and diversity of generated molecules.
# Usage notes
Each molecule should be stored in a separate SDF file. If all molecules are merged into a single SDF file, we also provide a Python script: split_sdf.py
### Required Configuration

| Parameter         | Description                                  | Example Value                                           |
|-------------------|----------------------------------------------|----------------------------------------------------------|
| `input_directory` | Directory containing input SDF files         | `"Diff-3D/data/QM9/sdf"`              |
| `output_directory`| Directory for saving output results          | `"Diff-3D/output/"`                 |
| `dataset`         | Dataset type (`QM9` or `Drug`)               | `"QM9"`                                                  |
| `forcefield`      | Force field type (`OPLS3` or `MMFF`)         | `"OPLS3"`                                                |
> 💡 Note: When using the OPLS3 force field, the ligprep_executable parameter must be specified to provide the path to the ligprep program.
### Reference File Paths

| Parameter           | Description                                                              | Example Value                      |
|---------------------|--------------------------------------------------------------------------|------------------------------------|
| train_smiles_file   | SMILES file of the training set used for novelty calculation             | "Diff-3D/data/data_reference/QM9/QM9_smiles.txt"          |
| training_set_file   | Training set file used for Fréchet Distance (FDC) calculation             | "Diff-3D/data/training_set_reference/QM9/smiles.smi"              |
| fpscores            | Fragment contribution score file used for SA score calculation           | "Diff-3D/data/fpscores.pkl.gz"         |

> 💡 **Note**  
There are two `train_smiles_file files` under `Diff-3D/data/data_reference/`, located in the `QM9/` and `Drug/` folders, respectively. These files contain SMILES strings of the training molecules and are used to calculate **novelty**.
>
>Similarly, the `training_set_file` is located under `Diff-3D/data/training_set_reference/`. Each dataset has its own file containing 10,000 SMILES strings randomly sampled from the training set, used for **FCD computation**.
### Evaluation Options

| Parameter                     | Description                                                             | Default Value |
|-------------------------------|-------------------------------------------------------------------------|---------------|
| calculate_qed                 | Whether to calculate QED (Quantitative Estimate of Drug-likeness)       | true          |
| calculate_sa                  | Whether to calculate SA (Synthetic Accessibility) score                 | true          |
| calculate_2d_properties       | Whether to calculate 2D structural properties                           | true          |
| calculate_stability           | Whether to calculate molecular stability                                | true          |
| calculate_fcd                 | Whether to calculate Fréchet ChemNet Distance and KL divergence scores  | true         |
| calculate_overall_diversity  | Whether to calculate overall molecular diversity                        | true          |
| calculate_scaffold_diversity | Whether to calculate scaffold diversity                                 | true          |
| calculate_geometry_metrics   | Whether to calculate 3D geometric metrics                                | true          |
# Example Usage

#### Step 1: Create and activate the Conda environment

There are two ways to set up the environment:

**Option 1 (Recommended):**
```bash
conda env create -f environment.yaml
conda activate diffmol
```
**Option 2 (Recommended):**
```bash
conda create -n diffmol python=3.8
conda activate diffmol
conda install -c conda-forge rdkit
conda install -c jrakhmonov pytorch
pip install guacamol
```
#### Step 2: Configure and run

Once your environment is ready, configure your `config.json` file with the required settings.

Then, simply run:

```bash
python main.py
```

