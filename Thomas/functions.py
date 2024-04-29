from rdkit import Chem
from rdkit.Chem import Draw
import time
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from bokeh.plotting import figure, show

def list_compiler ():
    df = pd.read_csv('Thomas/abundance.txt'
                    , sep='\t'
                    , header=None
                    , names=['Atom', 'Mass', 'Percentage'])
    mass = df['Mass'].tolist()
    mass = [float(m) for m in mass]
    abundance_percent = df['Percentage'].tolist()
    abundance_percent = [float(ap) for ap in abundance_percent]
    abundance = []
    for percent in abundance_percent:
        abundance.append(percent/100)
    isotopes = df['Atom'].tolist()
    return mass, abundance, isotopes

#Turn SMILEs representation into list of atomic symbols

mol_smi = input('Enter SMILEs: ')
mol_without_Hs = Chem.MolFromSmiles(mol_smi)
mol = Chem.AddHs(mol_without_Hs)


def atom_list_generator(mol):
    list_atoms_preionisation = []
    for atom in mol.GetAtoms():
        list_atoms_preionisation.append(atom.GetSymbol())
    return list_atoms_preionisation

def ionisation_method(list_atoms_preionisation):
    list_atoms = list_atoms_preionisation
    if 'H' in list_atoms:
        #Check that there is in fact a proton to remove
        list_atoms.remove('H')
    return list_atoms

print(ionisation_method(mol))