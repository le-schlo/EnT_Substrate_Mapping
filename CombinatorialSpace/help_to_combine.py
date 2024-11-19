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

def encode_mono_substitution_azoles(mol,
                                    patterns: List[str]
                                    ) -> List[int]:
    '''
    Given a MOL, return a list of 3 elements, where each element is 1 if the substitution is present, 0 otherwise

    order corresponds to carbons of thiazole: s1 c2 n3 c4 c5 [c2, c4, c5]

    :param mol: MOL object
    :param patterns: list of SMARTS patterns for the different substitution patterns
    :return: list of one hot encoded position of substitution
    '''

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
        substitution = [1,1,1]

    return substitution


def encode_mono_substitution_monocycles(mol,
                                    patterns: List[str]
                                    ) -> List[int]:
    '''
    Given a MOL, return a list of 3 elements, where each element is 1 if the substitution is present, 0 otherwise

    :param mol: MOL object
    :param patterns: list of SMARTS patterns for the different substitution patterns
    :return: list of one hot encoded position of substitution
    '''

    mono_2 = Chem.MolFromSmarts(patterns[0])
    mono_3 = Chem.MolFromSmarts(patterns[1])

    if mol.HasSubstructMatch(mono_2):
        substitution = [1, 0]
    elif mol.HasSubstructMatch(mono_3):
        substitution = [0, 1]

    else:
        substitution = [1, 1]

    return substitution


