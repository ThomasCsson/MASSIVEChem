from rdkit import Chem
from rdkit.Chem import Draw, AllChem

import pandas as pd

import base64

import os

from bokeh.plotting import figure, show, row
from bokeh.models.tickers import FixedTicker
from bokeh.layouts import row, column
from bokeh.io import show
from bokeh.models import ColumnDataSource, HTMLTemplateFormatter, WheelPanTool, WheelZoomTool, BoxAnnotation, CustomJS
from bokeh.models.widgets import DataTable, TableColumn

def empty_file_path(search_directory='.'):

     #---------------------------------------------------------------------------------------------#
    '''
    empty_file_path()
    
    Input: search_directory, which specifies to check in the current directory
    
    Output: - creates a new file called molecule_image.png to store an image later
            - returns the relative path to the file
    '''
    #---------------------------------------------------------------------------------------------#

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
    
    
    for root, dirs, files in os.walk(search_directory):

        #checks all the file names in the directory
        if filename in files:

            #return file path
            return os.path.join(root, filename)
    
    return None

def save_molecule_image_to_file(mol_smi, file_path, show_Hs=False, show_3D = False):

    #---------------------------------------------------------------------------------------------#
    '''
    save_molecule_image_to_file(mol_smi)
    
    Input: molecule under SMILEs representation
    
    Output: image of the molecule in the molecule_image.png file

    Functionnality: - if show_Hs= True:
                        shows all the hydrogens of the molecule
                    - if show_3D= True:
                        shows the molecule in 3D and the chirality
    '''
    #---------------------------------------------------------------------------------------------#


    # Generate the image from the molecule
    mol = Chem.MolFromSmiles(mol_smi)

    # Adds the hydrogens to the molecule if specified
    if show_Hs:
        mol = Chem.AddHs(mol)

    # Show the molecule in 3D if specified
    if show_3D:
        mol = AllChem.EmbedMolecule(mol)

    image = Draw.MolToImage(mol)

    # Save the image to a file
    image.save(file_path)

def mol_web_show(image_url):

    #---------------------------------------------------------------------------------------------#
    '''
    mol_web_show(image_url)
    
    Input: path of the molecule image
    
    Output: image of the molecule in bokeh

    '''
    #---------------------------------------------------------------------------------------------#


    # Creating a Bokeh figure to display the molecule
    p = figure(width=400, height=400,toolbar_location=None, x_range=(0, 1), y_range=(0, 1))
    p.image_url(url=[image_url], x=0, y=1, w=1, h=1)

    # Hide grid lines and axes
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    p.xaxis.visible = False
    p.yaxis.visible = False

    show(p)

    return p
