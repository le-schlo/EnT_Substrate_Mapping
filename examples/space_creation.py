from CombinatorialSpace.help_to_combine import get_functional_groups, create_space

path_to_fg_list = 'CombinatorialSpace/functional_groups.txt'
functional_groups = get_functional_groups(path_to_fg_list)

##Example Skipt to generate monosubstituted Thiophenes
smiles = 'C1=CC=CS1'
monosubstituted_thiophenes = create_space(functional_groups, smiles, 's1cccc1', num_iterations=1)

##Example Skipt to generate disubstituted Thiazoles
smiles = 'C1=CN=CS1'
disubstituted_thiazoles = create_space(functional_groups, smiles, 's1cnc1', num_iterations=2)

##Example Skipt to generate trisubstituted Furans
smiles = 'C1=CC=CO1'
trisubstituted_furans = create_space(functional_groups, smiles, 'o1ccc1', num_iterations=3)

##Example Skipt to generate mono- and disubstituted Oxazoles
smiles = 'C1=CC=CO1'
mono_and_disubstituted_oxazoles = create_space(functional_groups, smiles, 'o1ccc1', num_iterations=2, keep_all=True)

