def molecule_list_generator(mol) -> list[str]:
    #---------------------------------------------------------------------------------------------#
    '''
    molecule_list_generator(mol)
    
    Input: molecule under MOL representation
    
    Output: list containing the atomic symbol of each atom in the input molecule
    '''
    #---------------------------------------------------------------------------------------------#
    if not mol:
        raise ValueError('Enter a non-empty input')
    mol_test = mol.GetAtoms()
    if mol_test is None:
        raise ValueError('Invalid MOL representation')
    list_atoms = []
    for atom in mol.GetAtoms():
        list_atoms.append(atom.GetSymbol())
    return list_atoms

print(molecule_list_generator('CCXC'))