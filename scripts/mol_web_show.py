from bokeh.models import TextInput, Button, CustomJS, ColumnDataSource
from bokeh.plotting import figure, show, column
from rdkit import Chem
from rdkit.Chem import Draw
import os

def mol_web_show(mol_smi,show_Hs= False):
        # name of the file name to create
    filename = 'molecule_image.png'

    #finds the current directory
    current_directory = os.getcwd()

    #creates a file path to the current directory
    filepath = os.path.join(current_directory, filename)

    #checks if the file already exists
    if not os.path.exists(filepath):

        #if no, creates the path
        with open(filepath, 'a'):
            pass
    else:

        #else pass
        pass

    # Generate the image from the molecule
    mol = Chem.MolFromSmiles(mol_smi)
    # Adds the hydrogens to the molecule if specified
    if show_Hs:
        mol = Chem.AddHs(mol)

    image = Draw.MolToImage(mol)

    # Save the image to a file
    image.save(filepath)
    print(filepath)

    source = ColumnDataSource(data={'url': ['molecule_image.png']})

    # Create a text input for specifying the SMILES of the molecule
    smiles_input = TextInput(placeholder='Enter SMILES', title='SMILES')


     # Callback function to update the molecule image based on the input SMILES
    update_callback = CustomJS(args=dict(source=source, smiles_input=smiles_input), code="""
        var new_smiles = smiles_input.value;
        var xhr = new XMLHttpRequest();
        xhr.open("POST", "/update_image", true);
        xhr.setRequestHeader("Content-Type", "application/json");
        xhr.onreadystatechange = function () {
            if (xhr.readyState == 4 && xhr.status == 200) {
                // Upon successful response, update the image URL
                source.data['url'] = [xhr.responseText];
                source.change.emit();  // Trigger change event to update the plot
            }
        };
        xhr.send(JSON.stringify({smiles: new_smiles}));
    """)


    # Create a "Submit" button to trigger the update of the plot
    submit_button = Button(label='Submit', button_type='success')
    submit_button.js_on_click(update_callback)

    # Create a ColumnDataSource to hold the image URL
    source = ColumnDataSource(data={'url': ['molecule_image.png']})


    # Creating a Bokeh figure to display the molecule
    p = figure(width=350, height=350,toolbar_location=None, x_range=(0, 1), y_range=(0, 1))
    p.image_url(url=[filepath], x=0, y=1, w=1, h=1)

    # Hide grid lines and axes
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    p.xaxis.visible = False
    p.yaxis.visible = False

    layout = column(p, smiles_input, submit_button)
    
    return layout


mol_smi = input('MOL:  ')
show(mol_web_show(mol_smi))


import unittest

"""class TestMolWebShow(unittest.TestCase):
    def test_mol_web_show_without_Hs_3D(self):
        # Test the function without adding hydrogens
        p = mol_web_show('CCO')
        self.assertIsNotNone(p)

    def test_mol_web_show_with_Hs(self):
        # Test the function with adding hydrogens
        p = mol_web_show('CCO', show_Hs=True)
        self.assertIsNotNone(p)


if __name__ == '__main__':
    unittest.main()"""
