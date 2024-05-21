from rdkit import Chem

def atom_present_list(mol_smi):
    #---------------------------------------------------------------------------------------------#
    '''
    atom_present_list(mol_smi)
    
    Input: molecule under MOL representation
    
    Output: list of atoms present in the molecule with their respective count
    '''
    #---------------------------------------------------------------------------------------------#
    mol, list_atoms,list_atoms = Chem.MolFromSmiles(mol_smi),[],[]
    for atom in mol.GetAtoms():
        list_atoms.append(atom.GetSymbol())
    for atom in mol.GetAtoms():
        element = (f'{atom.GetSymbol()} : {list_atoms.count(atom.GetSymbol())}')   
        if element not in list_atoms:
            list_atoms.append(element)
    return list_atoms