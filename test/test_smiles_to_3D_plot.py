from rdkit import Chem
from rdkit.Chem import AllChem
import tempfile
from xyz2graph import MolGraph, to_plotly_figure
import panel as pn

def smiles_to_3D_plot(mol_smi):

     #---------------------------------------------------------------------------------------------#
    '''
    smiles_to_3D_plot(smiles)
    
    Input: molecule is a SMILEs format

    Output: panel graph of the molecule in 3D and interactive
    '''
    #---------------------------------------------------------------------------------------------#

    if not mol_smi:
        raise ValueError('Invalid SMILEs input')

    # Generate 3D coordinates from SMILES string
    mol2 = Chem.MolFromSmiles(mol_smi)

    if not mol2:
        raise ValueError('Invalide SMILEs input')
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

    return fig

import unittest

class TestSmilesTo3DPlot(unittest.TestCase):
    def test_valid_smiles(self):
        
        mol_smi = "CCO" 

        result = smiles_to_3D_plot(mol_smi)

        self.assertIsNotNone(result)

    def test_invalid_smiles(self):
        
        with self.assertRaises(ValueError):
            smiles_to_3D_plot("")  

    def test_nonexistent_smiles(self):
        
        with self.assertRaises(ValueError):
            smiles_to_3D_plot("invalid_smiles")    

    def test_plot_data(self):
        
        mol_smi = "CCO"  

        result = smiles_to_3D_plot(mol_smi)

        self.assertIn('data', result)  

    def test_plot_layout(self):
        
        mol_smi = "CCO"  

        result = smiles_to_3D_plot(mol_smi)

        self.assertIn('layout', result)  

    def test_plot_layout_attributes(self):
        
        mol_smi = "CCO" 

        result = smiles_to_3D_plot(mol_smi)
        layout = result['layout']
        
        self.assertIn('title', layout)  
        self.assertIn('scene', layout)
if __name__ == '__main__':
    unittest.main()


