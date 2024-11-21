# Accelerated Discovery in Dearomative Energy Transfer Catalysis through Data-Guided Reaction Screening

This repository contains the code for the creation and prediction of the virtual chemical space for the dearomative cycloaddition to heterocycles as described in this paper: TBA

<!---
<p align="center">
<img src="TOC.png" width="60%" />
</p>
-->

## Installation
For installation run
```
git clone --recurse-submodules https://github.com/le-schlo/EnT_Substrate_Mapping.git

pip install -r requirements.txt

#Additionally you need to go to the EasyChemML directory and install the necessary dependencies separately.
cd EnTdecker/EasyChemML
pip install ./
```

## Usage
The scripts to run the different parts of the workflow are located in the `examples` directory.
A typical workflow would involve the following steps:
1. Create the virtual chemical space: `examples/space_creation.py`
   - Define the core molecule (and matching SMARTS for decoration) for which the virtual library should be created. You can use the functional_groups in `CombinatorialSpace/functional_groups.txt` or define your own.
   - Save your virtual library as a `.csv` file to the `Data/virtual_libraries` directory. The naming convention should be `mol_substitution.csv`. Where `mol` is the name of the core molecule (e.g., thiophene) and `substitution` is the degree of substitution (e.g., monosubstituted).
   <br/><br/>
2. Obtain predictions: `examples/predict_et.py` and `examples/predict_sp.py`:
   - Download the pre-trained EnTdecker models from the [Zenodo](https://zenodo.org/records/10391170) repository and place them in the respective directories in the `Models` directory.
   - Alternatively you can also train / retrain your own models. For this please refer to the instructions in the [EnTdecker](https://github.com/le-schlo/EnTdecker) repository.
   - Run the scripts for prediction. The predictions should be save to `Data/et_predictions/predictions_mol_substitution.csv` for triplet energy prediction and `Data/sp_predictions/predictions_mol_substitution.csv` for spn population prediction.
   <br/><br/>
3. Create images for the interactive analysis (Optional): `examples/image_creation.py`
   - This script will create images. The user has two options controllable with the `type_of_image` variable
     - `structure`: The images will be the structures of the molecules generated with RDKit
     - `spin_population`: The images will be a heat map of predicted spin population of the molecules
     - The images should be saved to the `Data/images/idx.jpeg` Where `idx` is the index of the molecule in the virtual library. A dummy image named `default.png` is provided in the `Data/images` directory if no image for the molecule was created.
<br/><br/>
4. Run the interactive analysis: `examples/interactive_analysis.py`
   - Using the previously generated files an interactive 3D plot can be generated.
   - The representation of the molecules can be chosen with the `representation` variable:
     - `ECFP`: The molecules are represented by their ECFP fingerprints generated with RDKit
     - `MACCS`: The molecules are represented by their MACCS fingerprints generated with RDKit
     - `MACCS+ECFP`: The molecules are represented by the concatenation of their ECFP and MACCS fingerprints
   - The parameters for the dimensionality reduction with _UMAP_ can be adjusted:
     - `n_neighbors` default is 60. 
     - `min_dist` default is 0.15.
     - `metric` default is `euclidean`. 

## Citation
```
@article{TBA}
```
