from xyz2graph import MolGraph, to_networkx_graph, to_plotly_figure
from plotly.offline import offline

from rdkit import Chem
from rdkit.Chem import AllChem

import plotly.graph_objects as go
from bokeh.io import output_file, show
from bokeh.layouts import layout
from bokeh.models import Div
from plotly.io import to_html

def smiles_to_xyz(smiles, output_file):
    # Generate 3D coordinates from SMILES string
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)  # Add hydrogens for better geometry optimization
    AllChem.EmbedMolecule(mol, randomSeed=42)  # Embed the molecule in 3D space
    AllChem.MMFFOptimizeMolecule(mol)  # Optimize the geometry using MMFF94 force field

    # Write coordinates to XYZ file
    with open(output_file, 'w') as f:
        f.write(str(mol.GetNumAtoms()) + '\n\n')  # Write number of atoms
        for atom in mol.GetAtoms():
            pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
            f.write('{} {} {} {}\n'.format(atom.GetSymbol(), pos.x, pos.y, pos.z))

# Example usage:
smiles = "CCC(=O)OOCCC"  # Ethanol molecule
output_file = "ethanol.xyz"
smiles_to_xyz(smiles, output_file)


# Create the MolGraph object
mg = MolGraph()

# Read the data from the .xyz file
mg.read_xyz('ethanol.xyz')

# Create the Plotly figure object
fig = to_plotly_figure(mg)

# Plot the figure
#offline.plot(fig)


import numpy as np
import plotly.graph_objs as go

import panel as pn

pn.extension("plotly")

xx = np.linspace(-3.5, 3.5, 100)
yy = np.linspace(-3.5, 3.5, 100)
x, y = np.meshgrid(xx, yy)
z = np.exp(-((x - 1) ** 2) - y**2) - (x**3 + y**4 - x / 5) * np.exp(-(x**2 + y**2))

surface=go.Surface(z=z)
fig = go.Figure(data=[surface])

fig.update_layout(
    title="Plotly 3D Plot",
    width=500,
    height=500,
    margin=dict(t=50, b=50, r=50, l=50),
)

plotly_pane = pn.pane.Plotly(fig)
plotly_pane