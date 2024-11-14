from rdkit import Chem
from rdkit.Chem import AllChem
import tempfile
from xyz2graph import MolGraph

def smiles_to_3D_plot(mol_smi):

     #---------------------------------------------------------------------------------------------#
    '''
    smiles_to_3D_plot(smiles)
    
    Input: molecule is a SMILEs format

    Output: panel graph of the molecule in 3D and interactive
    '''
    #---------------------------------------------------------------------------------------------#

    if not mol_smi:
        raise ValueError("Enter a non-empty string as argument")

    # Generate 3D coordinates from SMILES string
    mol2 = Chem.MolFromSmiles(mol_smi)

    if not mol2:
        raise ValueError('Enter a correct SMILEs')
    
    mol2 = Chem.AddHs(mol2)  # Add hydrogens for better geometry optimization
    AllChem.EmbedMolecule(mol2, randomSeed=42)  # Embed the molecule in 3D space
    AllChem.MMFFOptimizeMolecule(mol2)  # Optimize the geometry using MMFF94 force field

    # Create a temporary file to store the 
    with tempfile.NamedTemporaryFile(delete=False) as tmp:
        tmp.write(f"{mol2.GetNumAtoms()}\n\n".encode('utf-8'))  # Write number of atoms
        for atom in mol2.GetAtoms():
            pos = mol2.GetConformer().GetAtomPosition(atom.GetIdx())
            tmp.write(f"{atom.GetSymbol()} {pos.x} {pos.y} {pos.z}\n".encode('utf-8'))
        tmp_path = tmp.name

    # Create the MolGraph object
    mg = MolGraph()
    
    # Read the data from the temporary file
    mg.read_xyz(tmp_path)

    # Create the Plotly figure object
    fig = mg.to_plotly()

    return fig






