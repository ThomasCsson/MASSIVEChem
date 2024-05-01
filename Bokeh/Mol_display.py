
from bokeh.plotting import figure, show
from rdkit.Chem import Draw, AllChem
from rdkit import Chem
from bokeh.plotting import row
import numpy as np
import matplotlib as plt


input_mol = input('SMILES: ')

# Function to generate RDKit molecule image and save to file
def save_molecule_image_to_file(smi, file_path, show_Hs=False, show_3D = False):
    # Generate the image from the molecule
    mol = Chem.MolFromSmiles(smi)

    # Adds the hydrogens to the molecule if specified
    if show_Hs:
        mol = Chem.AddHs(mol)

    # Show the molecule in 3D if specified
    if show_3D:
        mol_3D = AllChem.EmbedMolecule(mol)

    image = Draw.MolToImage(mol)

    # Save the image to a file
    image.save(file_path)

output_file_path = "molecule_image.png"
save_molecule_image_to_file(input_mol, output_file_path, False, False)

def mol_web_show(image_url):

    # Creating a Bokeh figure
    p = figure(width=400, height=400,toolbar_location=None, x_range=(0, 1), y_range=(0, 1))
    p.image_url(url=[image_url], x=0, y=1, w=1, h=1)

    # Hide grid lines and axes
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    p.xaxis.visible = False
    p.yaxis.visible = False

    return p



image_url = "/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/molecule_image.png"

#show(mol_web_show(image_url))

def final_layout(p1, p2):
    layout = row(p1,p2)
    show(layout)

input_mol = input('SMILEs : ')

functional_groups_smiles = {
    'Alcohol': 'CO',
    'Aldehyde': 'CC=O',
    'Ketone': 'CC(=O)C',
    'Carboxylic Acid': 'CC(=O)O',
    'Ester': 'CC(=O)OC',
    'Ether': 'CCOCC',
    'Amide': 'CC(=O)N',
    'Amine': 'CN',
    'Nitrile': 'CC#N',
    'Chloride': 'CCl',
    'Bromide': 'CBr',
    'Fluoride': 'CF',
    'Iodide': 'CI',
    'Alkene': 'C=C',
    'Alkyne': 'C#C',
    'Imine': 'C=NC',
    'Amino acid': 'CC(N)C(=O)O',
    'Thiol': 'CS',
    'Epoxyde': '[*]C1([*])OC1([*])[*]',
    'Sulfides': 'CSC',
    'Acyl Chloride': 'CC(=O)Cl',
    'Anhydride': 'CC(=O)OC(=O)C',
    'Nitro': 'C[N+](=O)[O-]',
    'Enamine': 'C=CN',
    'Imide': 'C(=O)NC(=O)C',
    'Azide': 'CNNN',
    'Enol': 'C=C(O)C',
    'Hemiacetal': 'CC(O)C',
    'Carbonate': 'OC(=O)O',
    'Disulfide': 'CSSC',
    'Sulfoxide': 'CS(=O)C',
    'Sulfone': 'CS(=O)(=O)C',
    'Sulfonic acid': 'CS(=O)(=O)O',
    'Sulfanote ester': 'CS(=O)(=O)OC',
    'Thioester': 'C(=O)SC',
    'Phosphine': 'CP',
    'Phosphate ester': 'COP(=O)(O)O'
}

def subgroup_value(mol_smi, dict_functional_groups):
    list_contained_subgroups, list_contained_subgroups_values = [], []
    mol = Chem.MolFromSmiles(mol_smi)

    for SMILES, value in dict_functional_groups.items():
        substruct = Chem.MolFromSmiles(dict_functional_groups[SMILES])  # Get the SMILES string from the dictionary

        if mol.HasSubstructMatch(substruct):
            list_contained_subgroups.append(SMILES)
            list_contained_subgroups_values.append(value)
            indices = mol.GetSubstructMatch(substruct)
            print(indices)
        else:
            list_contained_subgroups.append(f'NOT {SMILES}')

    if 'Aldehyde' and 'Ketone' in list_contained_subgroups:
        mol_aldehyde = Chem.MolFromSmiles('CC=O')
        mol_ketone = Chem.MolFromSmiles('C(=O)C')

        list_aldehyde_index = mol.GetSubstructMatch(mol_aldehyde)
        list_ketone_index = mol.GetSubstructMatch(mol_ketone)

        for x in list_aldehyde_index:
            if x in list_ketone_index:
                list_contained_subgroups[1] = 'NOT aldehyde'

    if 'Carboxylic acid' and 'Aldehyde' and 'Alcohol' in list_contained_subgroups:

        mol_acid = Chem.MolFromSmiles('CC(=O)O')
        mol_aldehyde = Chem.MolFromSmiles('CC=O')
        mol_alcohol = Chem.MolFromSmiles('CO')

        list_aldehyde_index = mol.GetSubstructMatch(mol_aldehyde)
        list_acid_index = mol.GetSubstructMatch(mol_acid)
        list_alcohol_index = mol.GetSubstructMatch(mol_alcohol)

        for x in list_acid_index:
            if x in list_aldehyde_index:
                list_contained_subgroups[1] = 'NOT Aldehyde'
            if x in list_alcohol_index:
                list_contained_subgroups[0] = 'NOT Alcohol'

    if 'Carboxylic acid' and 'Ester' and 'Ether' in list_contained_subgroups:

        mol_acid = Chem.MolFromSmiles('CC(=O)O')
        mol_ester = Chem.MolFromSmiles('CC(=O)OC')
        mol_ether = Chem.MolFromSmiles('COC')

        list_acid_index = mol.GetSubstructMatch(mol_acid)
        list_ester_index = mol.GetSubstructMatch(mol_ester)
        list_ether_index = mol.GetSubstructMatch(mol_ether)

        for x in list_ester_index:
            if x in list_acid_index:
                list_contained_subgroups[3] = 'NOT Carboxylic Acid'
            if x in list_ether_index:
                list_contained_subgroups[5] = 'NOT Ether'

    if 'Ether' and 'Alcohol' in list_contained_subgroups:
        mol_alcohol = Chem.MolFromSmiles('CO')
        mol_ether = Chem.MolFromSmiles('COC')

        list_alcohol_index = mol.GetSubstructMatch(mol_alcohol)
        list_ether_index = mol.GetSubstructMatch(mol_ether)


    return list_contained_subgroups, list_contained_subgroups_values


print(subgroup_value(input_mol,functional_groups_smiles))


#alcohol ok
#aldehyde ok
#ketone ok