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
    
    if not mol_smi:
        return None
    
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

import unittest


import unittest
from bokeh.plotting import figure


class TestSmilesTo3DPlot(unittest.TestCase):

    def create_random_bokeh_figure(self):
        # Create a random Bokeh figure
        p = figure(plot_width=300, plot_height=300)
        p.circle([1, 2, 3, 4, 5], [6, 7, 2, 4, 5], size=10, color="navy", alpha=0.5)
        return p

    def test_smiles_to_3D_plot_output_type(self):
        # Test if the output type is correct
        mol_smi = "CCO"  # Example SMILES string
        last_plot = self.create_random_bokeh_figure()  # Random Bokeh figure
        double_graph = self.create_random_bokeh_figure()  # Random Bokeh figure
        result = smiles_to_3D_plot(mol_smi, last_plot, double_graph)
        self.assertIsInstance(result, pn.Row)  # Check if the output is a Row element

    def test_smiles_to_3D_plot_input_smiles(self):
        # Test with an invalid SMILES string
        mol_smi = "invalid_smiles"  # Example invalid SMILES string
        last_plot = self.create_random_bokeh_figure()  # Random Bokeh figure
        double_graph = self.create_random_bokeh_figure()  # Random Bokeh figure
        result = smiles_to_3D_plot(mol_smi, last_plot, double_graph)
        # Check if the output is not None (meaning the function didn't crash with invalid input)
        self.assertIsNotNone(result)

    def test_smiles_to_3D_plot_input_last_plot(self):
        # Test with None as the last plot
        mol_smi = "CCO"  # Example SMILES string
        last_plot = None  # None as the last plot
        double_graph = self.create_random_bokeh_figure()  # Random Bokeh figure
        result = smiles_to_3D_plot(mol_smi, last_plot, double_graph)
        # Check if the output is not None (meaning the function didn't crash with None as last plot)
        self.assertIsNotNone(result)

    def test_smiles_to_3D_plot_input_double_graph(self):
        # Test with None as the double graph
        mol_smi = "CCO"  # Example SMILES string
        last_plot = self.create_random_bokeh_figure()  # Random Bokeh figure
        double_graph = None  # None as the double graph
        result = smiles_to_3D_plot(mol_smi, last_plot, double_graph)
        # Check if the output is not None (meaning the function didn't crash with None as double graph)
        self.assertIsNotNone(result)

if __name__ == '__main__':
    unittest.main()






