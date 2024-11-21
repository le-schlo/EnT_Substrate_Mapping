import sys
import pandas as pd

home_dir = r'/path/to/project_directory/' #Change this to the directory your project is located at.
sys.path.append(home_dir) #Add the directory to your pythonpath, so that the functions from Combinatorial.help_to_combine can be imported.

from CombinatorialSpace.help_to_combine import get_functional_groups, create_space

path_to_fg_list = home_dir + 'CombinatorialSpace/functional_groups.txt'
functional_groups = get_functional_groups(path_to_fg_list)

##Example Skipt to generate monosubstituted Thiophenes
smiles = 'C1=CC=CS1'
monosubstituted_thiophenes = create_space(functional_groups, smiles, 's1cccc1', num_iterations=1)

# ##Example Skipt to generate disubstituted Thiazoles
# smiles = 'C1=CN=CS1'
# disubstituted_thiazoles = create_space(functional_groups, smiles, 's1cnc1', num_iterations=2)
#
# ##Example Skipt to generate trisubstituted Furans
# smiles = 'C1=CC=CO1'
# trisubstituted_furans = create_space(functional_groups, smiles, 'o1ccc1', num_iterations=3)
#
# ##Example Skipt to generate mono- and disubstituted Oxazoles
# smiles = 'C1=CC=CO1'
# mono_and_disubstituted_oxazoles = create_space(functional_groups, smiles, 'o1ccc1', num_iterations=2, keep_all=True)

df = pd.DataFrame(data={'smiles': monosubstituted_thiophenes})
df.to_csv(home_dir + 'Data/virtual_libraries/thiophene_monosubstituted.csv', index=False)

