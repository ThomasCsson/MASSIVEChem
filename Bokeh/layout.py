from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

def save_molecule_image_to_file(mol_smi, file_path, show_Hs=False, show_3D = False):

    #---------------------------------------------------------------------------------------------#
    '''
    save_molecule_image_to_file(mol_smi)
    
    Input: molecule under SMILEs representation
    
    Output: image of the molecule in the molecule_image.png file

    Functionnality: - if show_Hs= True:
                        shows all the hydrogens of the molecule
                    - if show_3D= True:
                        shows the molecule in 3D and the chirality
    '''
    #---------------------------------------------------------------------------------------------#


    # Generate the image from the molecule
    mol = Chem.MolFromSmiles(mol_smi)

    # Adds the hydrogens to the molecule if specified
    if show_Hs:
        mol = Chem.AddHs(mol)

    # Show the molecule in 3D if specified
    if show_3D:
        mol = AllChem.EmbedMolecule(mol)

    image = Draw.MolToImage(mol)

    # Save the image to a file
    image.save(file_path)

filepath = './molecule_image_5.png'
mol = 'CCC(=O)OCC(=O)OO'

save_molecule_image_to_file(mol, filepath)


