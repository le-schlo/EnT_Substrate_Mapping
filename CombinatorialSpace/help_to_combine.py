from rdkit import Chem
from typing import List

def get_functional_groups(path_to_fg: str
                          ) -> List[str]:
    '''
    Given a path to a file containing functional group SMILES separated by \n, return a list FGs
    :param path_to_fg: Path to file containing functional group SMILES separated by \n
    :return: List of functional group SMILES
    '''
    fg_smiles = []
    with open(path_to_fg, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            fg_smiles.append(line)

    return fg_smiles

def combine_with_heterocycle(fg: str,
                             heterocycle_smi: str,
                             heterocycle_smarts: str
                             ) -> List[str]:
    '''
    Provide a functional group SMILES and return the combined molecule with a heterocycle provided in SMILES format (e.g., 'C1=CN=CO1' for oxazole).
    Provide the SMARTS pattern of the heterocycle. To all carbon atoms in this heterocycle that are not substituted yet, the functional group will be attached.

    :param fg: SMILES of functional group
    :param heterocycle_smi: SMILES of heterocycle
    :param heterocycle_smarts: SMARTS of heterocycle
    :return: list of combines molecules in SMILES format
    '''


    heterocycle_mol = Chem.MolFromSmiles(heterocycle_smi)
    heterocycle_patt_mol = Chem.MolFromSmiles(heterocycle_smarts)
    fg_mol = Chem.MolFromSmiles(fg)

    # combine the heterocycle with the functional group and get the matched atoms of the heterocycle
    combined_mol = Chem.CombineMols(heterocycle_mol,fg_mol)
    hit_atoms = list(combined_mol.GetSubstructMatch(heterocycle_patt_mol))
    mult_fg = list(combined_mol.GetSubstructMatches(fg_mol))

    #Add the functional groups to the carbons of the heterocycle that are not substituted yet and create a list of all possible combinations
    list_of_combinations = []
    occupied_atoms = []
    for atom_idx in hit_atoms:
        if combined_mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 6 and combined_mol.GetAtomWithIdx(atom_idx).GetTotalNumHs() != 0 and atom_idx not in occupied_atoms:

            attachement_atom_1 = atom_idx
            attachement_atom_2 = mult_fg[-1][0]
            em = Chem.EditableMol(combined_mol)
            em.AddBond(attachement_atom_1,attachement_atom_2,Chem.BondType.SINGLE)
            nm = em.GetMol()
            Chem.SanitizeMol(nm)
            list_of_combinations.append(Chem.MolToSmiles(nm))
            occupied_atoms.append(attachement_atom_1)

    return list_of_combinations

def create_space(fg: List[str],
                 heterocycle_smi: str,
                 heterocycle_smarts: str,
                 num_iterations: int = 1,
                 keep_all: bool = False
                 ) -> List[str]:
    '''
    Create a combinatorial space by combining a list of functional groups with a heterocycle

    :param fg: List of SMILES of functional groups
    :param heterocycle_smi: SMILES of heterocycle
    :param heterocycle_smarts: SMARTS of heterocycle
    :param num_iterations: Number of iterations to combine the functional groups with the heterocycle. Corresponds to mono-, di-, tri- substituted molecules
    :param keep_all: If True, all generated molecules (e.g., monosubstituted and disusbstituted molecules) are kept, if False, only molecules o last iteration (e.g., disubstituted) are kept
    :return: List of combined molecules in SMILES format
    '''

    generated_molecules = []
    for i in fg:
        generated_molecules.extend(combine_with_heterocycle(i, heterocycle_smi, heterocycle_smarts))

    unique_molecules = list(set(generated_molecules))

    for i in range(num_iterations - 1):
        if not keep_all:
            generated_molecules = []
        for smi in unique_molecules:
            for i in fg:
                generated_molecules.extend(
                    combine_with_heterocycle(i, heterocycle_smi=smi, heterocycle_smarts=heterocycle_smarts))
        unique_molecules = list(set(generated_molecules))

    return unique_molecules
def encode_mono_substitution_azoles(smi: str,
                                    patterns: List[str] = ['n1[cH][cH][s,o,nH]c1',
                                                           'n1c[cH][s,o,nH][cH]1',
                                                           'n1[cH]c[s,o,nH][cH]1']
                                    ) -> List[int]:
    '''
    Given a SMILES, return a list of 3 elements, where each element is 1 if the substitution is present, 0 otherwise

    order corresponds to carbons of thiazole: s1 c2 n3 c4 c5 [c2, c4, c5]

    :param str: SMILES of molecule
    :param patterns: list of SMARTS patterns for the different substitution patterns
    :return: list of one hot encoded position of substitution
    '''
    mol = Chem.MolFromSmiles(smi)
    mono_2 = Chem.MolFromSmarts(patterns[0])
    mono_4 = Chem.MolFromSmarts(patterns[1])
    mono_5 = Chem.MolFromSmarts(patterns[2])

    if mol.HasSubstructMatch(mono_2):
        substitution = [1,0,0]
    elif mol.HasSubstructMatch(mono_4):
        substitution = [0,1,0]
    elif mol.HasSubstructMatch(mono_5):
        substitution = [0,0,1]

    else:
        raise ValueError('No substitution found. Only monosubstituted molecules are supported')

    return substitution


def encode_mono_substitution_monoheterocycles(smi: str,
                                    patterns: List[str] = ['[o,s,n]1c[cH][cH][cH]1',
                                                           '[o,s,n]1[cH]c[cH][cH]1']
                                    ) -> List[int]:
    '''
    Given a SMILES, return a list of 3 elements, where each element is 1 if the substitution is present, 0 otherwise

    order corresponds to carbons of heterocycle, e.g., thiophene: s1 c2 c3 c4 c5 [c2, c3]

    :param smi: Smiles of molecule
    :param patterns: list of SMARTS patterns for the different substitution patterns
    :return: list of one hot encoded position of substitution
    '''
    mol = Chem.MolFromSmiles(smi)
    mono_2 = Chem.MolFromSmarts(patterns[0])
    mono_3 = Chem.MolFromSmarts(patterns[1])

    if mol.HasSubstructMatch(mono_2):
        substitution = [1, 0]
    elif mol.HasSubstructMatch(mono_3):
        substitution = [0, 1]
    else:
        raise ValueError('No substitution found. Only monosubstituted molecules are supported')
    return substitution


