import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import SimilarityMaps
from rdkit.Chem import Draw



home_dir = '/path/to/project_directory/'
ring_class = 'thiophene'
substitution = 'monosubstituted'
type_of_image = 'structure' # Options are structure or 'spin_population'
                            # Caveat of spin_population is that the memory for the generation for the images can overflow and crash the script.
                            # In this case the script needs to be run in smaller batches.

data_path = home_dir + fr'\Data\virtual_libraries\{ring_class}_{substitution}.csv'
picFolder_path = home_dir + rf'\Data\images'
#Only required for spin_density images:
sp_path = home_dir + fr'\Data\sp_predictions\predictions_{ring_class}_{substitution}.csv'

df = pd.read_csv(data_path)

if not os.path.exists(picFolder_path):
    os.mkdir(picFolder_path)

if type_of_image == 'structure':
    for smi in df.smiles.to_list():
        d = Draw.rdMolDraw2D.MolDraw2DCairo(512, 512)
        d.drawOptions().rotate = 0
        d.drawOptions().bondLineWidth = 1
        d.DrawMolecule(Chem.MolFromSmiles(smi))
        d.FinishDrawing()
        idx = df[df['smiles'] == smi].index.to_numpy()[0]
        d.WriteDrawingText(f"{picFolder_path}/{idx}.jpeg")

elif type_of_image == 'spin_density':
    sp_data = pd.read_csv(sp_path)
    for i in range(len(sp_data)):
        smi = sp_data.True_smiles.to_list()[i]
        sp_list = sp_data.DensityArray.to_list()[i]
        string_list = sp_list.strip('[]')
        string_elements = string_list.split()
        sp_list = [int(element) for element in string_elements]
        m = Chem.MolFromSmiles(smi)
        densities = np.array(sp_list)
        count = df[df['smiles'] == smi].index.to_numpy()[0]
        fig = SimilarityMaps.GetSimilarityMapFromWeights(m,
                                                         densities,
                                                         colorMap='coolwarm',
                                                         contourLines=10,
                                                         colors='grey')
        fig.savefig(f"{picFolder_path}/{count}.jpeg",
                    bbox_inches='tight',
                    dpi=600)
        fig.clear()
