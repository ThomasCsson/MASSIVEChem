from rdkit import Chem
            
def atom_present_list(mol_smi):
    mol, list_atoms,list_out_final = Chem.MolFromSmiles(mol_smi),[],[]
    for atom in mol.GetAtoms():
        list_atoms.append(atom.GetSymbol())
    for atom in mol.GetAtoms():
        element = (f'{atom.GetSymbol()} : {list_atoms.count(atom.GetSymbol())}')   
        if element not in list_out_final:
            list_out_final.append(element)
    return list_out_final

print(atom_present_list('CCOCCOCCOCCOCC'))