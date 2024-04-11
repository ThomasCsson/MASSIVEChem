from rdkit import Chem

bonds = []
mol = Chem.MolFromSmiles('C(=O)CC#N')
for bond in mol.GetBonds():
    bonds.append(bond.GetBondType())
    print(bond.GetBondType())
print(bonds)
    