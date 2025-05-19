import pandas as pd
from rdkit import Chem

#--------------Define which data shall be visualised and set the corresponding paths-------------------#

ring_class = 'thiophene'
substitution = 'monosubstituted'

core_smiles = 'c1cccs1'
threshold = 62

home_dir = '/path/to/project_directory/'

path_to_et_datafile = home_dir + fr'Data/et_predictions/predictions_{ring_class}_{substitution}.csv'
path_to_sp_datafile = home_dir + fr'Data/sp_predictions/predictions_{ring_class}_{substitution}.csv'
path_to_save_candidates = home_dir + f'Data/potential_candidates_{ring_class}_{substitution}.csv'

#--------------Load the data, and check if the maximum spin density is located on the core-------------------#
data_et = pd.read_csv(path_to_et_datafile)
data_sp = pd.read_csv(path_to_sp_datafile)

data_sp = data_sp.rename(columns={'True_smiles': 'smiles'})
data = data_sp.merge(data_et[['smiles', 'e_t']], on='smiles', how='left')

data['DensityArray'] = data['DensityArray'].apply(lambda x: eval(x.replace(' ', ',')))

density_located_on_core = []

for index, row in data.iterrows():
    predicted_smiles = row['Predicted_smiles']
    density_array = data[data['smiles'] == predicted_smiles]['DensityArray'].values[0]

    pattern_to_match = Chem.MolFromSmiles(core_smiles)
    atom_index_of_core = list(Chem.MolFromSmiles(predicted_smiles).GetSubstructMatch(pattern_to_match))

    density_on_core_list = []
    density_on_substituent_list = []
    for i in range(len(density_array)):
        if i in atom_index_of_core:
            density_on_core_list.append(density_array[i])
        else:
            density_on_substituent_list.append(density_array[i])

    if max(density_on_core_list) >= max(density_on_substituent_list):
        density_on_core=1
    else:
        density_on_core=0

    density_located_on_core.append(density_on_core)

data['Density_on_core'] = density_located_on_core


#--------------Select the molecules that have a predicted triplet energy below the threshold and where the spin density is located on the ring-------------#

interesting_data = data[data['e_t'] < threshold]
interesting_data = interesting_data[interesting_data['Density_on_core'] == 1]
interesting_data.to_csv(path_to_save_candidates)
