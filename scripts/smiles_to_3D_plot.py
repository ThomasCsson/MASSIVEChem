from rdkit import Chem
from rdkit.Chem import AllChem
import tempfile
from xyz2graph import MolGraph, to_plotly_figure
import panel as pn

def smiles_to_3D_plot(mol_smi,last,double_graph):

     #---------------------------------------------------------------------------------------------#
    '''
    smiles_to_3D_plot(smiles)
    
    Input: molecule is a SMILEs format

    Output: panel graph of the molecule in 3D and interactive
    '''
    #---------------------------------------------------------------------------------------------#

    # Generate 3D coordinates from SMILES string
    mol2 = Chem.MolFromSmiles(mol_smi)
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
    fig = to_plotly_figure(mg)

    plotly_pane = pn.pane.Plotly(fig)

    p1_pane = pn.pane.Bokeh(last)
    p2_pane = pn.pane.Bokeh(double_graph)

    layout12 = pn.Column( plotly_pane, p1_pane)
    final = pn.Row(p2_pane, layout12)

    return final

# Example usage
input_mol = input('MOL:  ')
smiles_to_3D_plot(input_mol).show()




