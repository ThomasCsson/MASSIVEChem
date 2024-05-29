from rdkit import Chem
from rdkit.Chem import Draw

from io import BytesIO
import base64
from bokeh.models import Div


def mol_web_show(mol_smi):

    #---------------------------------------------------------------------------------------------#
    '''
    mol_web_show(mol_smi)
    
    Input: SMILEs of a molecule
    
    Output: image of the molecule as a bokeh plot
    '''
    #---------------------------------------------------------------------------------------------#

    if not mol_smi:
        raise ValueError('Incorrect SMILEs input')
    
    # Generate the image from the molecule
    mol = Chem.MolFromSmiles(mol_smi)

    if mol is None:
        raise ValueError('Incorrect SMILEs input')

    #Draws the image
    image = Draw.MolToImage(mol)

    #stocks the image in a base64 format
    buffered = BytesIO()
    image.save(buffered, format="PNG")
    image_base64 = base64.b64encode(buffered.getvalue()).decode("utf-8")
    image_url = f"data:image/png;base64,{image_base64}"

    # Create a Div element to display the image
    img_div = Div(text=f'<img src="{image_url}" style="width:350px;height:350px;">')
   
    return img_div
    
