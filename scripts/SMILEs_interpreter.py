from rdkit import Chem

def SMILEs_interpreter(mol_smi):
    #---------------------------------------------------------------------------------------------#
    '''
    SMILEs_interpreter(mol_smi)

    Input: molecule under SMILEs representation
    
    Output: molecule under MOL representation
    '''
    #---------------------------------------------------------------------------------------------#

    #Checks that the SMILEs input is correct/allowed

    mol_without_Hs = Chem.MolFromSmiles(mol_smi)

    if mol_without_Hs is None:
        print('')
        print("Invalid SMILEs input.")
        print('Please try again with a different SMILEs.')
        exit()

    mol = Chem.AddHs(mol_without_Hs)

    return mol