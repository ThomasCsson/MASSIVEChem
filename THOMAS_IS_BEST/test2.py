from rdkit import Chem
from rdkit.Chem import Draw
input = Chem.MolFromSmiles('CSC')



pat = Chem.MolFromSmarts("C[Sh1]")


print(input.HasSubstructMatch(pat))


    
        