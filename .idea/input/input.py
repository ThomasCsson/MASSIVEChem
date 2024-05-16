
from rdkit import Chem
from bokeh.models import Button, Div, CustomJS
from bokeh.layouts import column, row
from bokeh.plotting import show


def calculate_insaturation(mol_smile) -> int:
    #---------------------------------------------------------------------------------------------#
    '''
    This function takes as input a string representing the SMILES notation of a molecule. 
    Using the RDKit library, it creates a representation of the molecule, adds hydrogens and then calculates the unsaturation level of the molecule. 
    The level of unsaturation is determined by counting the number of carbon, nitrogen and halogen atoms, then applying a specific formula. 
    Finally, the function displays the unsaturation level and the image of the molecule, and indicates the time required to perform the function.
    '''
    #---------------------------------------------------------------------------------------------#

    mol_1 = Chem.MolFromSmiles(mol_smile)
    mol = Chem.AddHs(mol_1)

    C, N, HX = 0, 0, 0
    halogens_hydrogen = ['F', 'Cl', 'Br', 'I', 'At', 'H']
    for atom in mol.GetAtoms():
        atom_sym = atom.GetSymbol()
        if atom_sym == 'C':
            C += 1
        elif atom_sym == 'N':
            N += 1
        elif atom_sym in halogens_hydrogen:
            HX += 1

    insaturation = C + 1 + (N - HX) / 2

    return insaturation


def customJS():
    # Callback function to update the displayed information
    update_callback = CustomJS(code="""
    // Call the Python function to retrieve custom information
    var custom_info = Bokeh.index[model_id].get_model_by_name('custom_info');
    custom_info.text = "Unsaturation: " + calculate_insaturation(mol_smi);
    """)

    # Create a button
    button = Button(label="Update Information", button_type="success")
    button.js_on_click(update_callback)

    # Create a Div element to display the custom information
    info_div = Div(text="Custom information will be displayed here.", name="unsaturation")

    # Create a layout containing the button and the Div element
    layout = column(row(button), row(info_div))

    # Show the layout
    return layout


mol_smi = 'CCCO'

show(customJS())




