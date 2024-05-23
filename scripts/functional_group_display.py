from bokeh.models import ColumnDataSource, HTMLTemplateFormatter
from bokeh.models.widgets import DataTable, TableColumn
from bokeh.plotting import show
import base64
from io import BytesIO
from bokeh.models import ColumnDataSource, HTMLTemplateFormatter
from bokeh.models.widgets import DataTable, TableColumn
from bokeh.plotting import show
from rdkit import Chem
from rdkit.Chem import Draw


def functional_group_display(contained_functional_groups):

    #---------------------------------------------------------------------------------------------#
    '''
    functional_group_display(contained_functional_groups)
    
    Input: list of all contained functional groups in the molecule

    Output: bokeh table with the contained functional groups and their image
    '''
    #---------------------------------------------------------------------------------------------#

    #dictionnary with functional groups and associated SMARTs
    functional_groups_smarts = {
        'Alcohol': 'C[Oh1+0]',
        'Aldehyde': 'C[Ch1]=O',
        'Ketone': 'CC(=O)C',
        'Carboxylic Acid': 'CC(=O)[Oh1]',
        'Ester': 'CC(=O)[Oh0]',
        'Ether': '*[Oh0]*',
        'Amide': 'C(=O)N',
        'Amine': '[C][N]',
        'Nitrile': 'C#N',
        'Chloride': 'Cl',
        'Bromide': 'Br',
        'Fluoride': 'F',
        'Iodide': 'I',
        'Alkene': 'C=C',
        'Alkyne': 'C#C',
        'Imine': 'C=N*',
        'Amino acid': '[Nh2][Ch1*]C(=O)O',
        'Proline': '[Nh1][Ch1*]C(=O)O',
        'Thiol': '[Sh1]',
        'Sulfide': '*[Sh0]*',
        'Acyl Chloride': 'CC(=O)Cl',
        'Anhydride': '*[Ch0](=O)O[Ch0](=O)*',
        'Nitro': 'C[N+](=O)[O-]',
        'Enamine': 'C=C[Nh0]',
        'Enamine2': 'C=C[Nh1]',
        'Enamine3': 'C=C[Nh2]',
        'Imide': 'C(=O)NC(=O)*',
        'Azide': 'CNNN',
        'Enol': 'C=C([Oh1])C',
        'Hemiacetal': 'CC(O)(O)C',
        'Carbonate': '[Oh0]C(=O)[Oh0]',
        'Carbonate2': '[Oh1]C(=O)[Oh1]',
        'Disulfide': 'CSSC',
        'Sulfoxide': 'CS(=O)C',
        'Sulfone': '*[So2](=O)(=O)*',
        'Sulfonic acid': '*S(=O)(=O)[Oh1]',
        'Thioester': 'C(=O)S*',
        'Phosphine': '*[Po0](*)*',
        'Phosphate': '*OP(=O)(O)O',
        'Benzene': 'c1ccccc1',
        'Peroxide':'C[Oh0][Oh0]C'
    }

    #initiate empty variables
    present_group_smarts = []    
    present_group_images_base64 = []
    
    #appends the smarts of the contained functional groups
    for i,j in functional_groups_smarts.items():
        for x in contained_functional_groups:
            if x == i:
                present_group_smarts.append(j)

    #converts the smarts to images in base64 format
    for x in present_group_smarts:

        #converts SMARTs to SMILEs for the images to be nicer
        mol_x = Chem.MolFromSmarts(x)
        mol_smi = Chem.MolToSmiles(mol_x)
        mol = Chem.MolFromSmiles(mol_smi)

        if mol:

            #creates the image
            image = Draw.MolToImage(mol)

            # Convert the image to base64 format
            buffered = BytesIO()
            image.save(buffered, format="PNG")
            image_base64 = base64.b64encode(buffered.getvalue()).decode("utf-8")

            #appends the image to an empty dictionnary
            present_group_images_base64.append(image_base64)

    #creates a dictionnary that links the name of the functional group to its image
    data = dict(
        groups=contained_functional_groups,
        images=present_group_images_base64
    )
    source = ColumnDataSource(data)

    #template for the bokeh table that read the base64 format 
    template = """
    <div>
        <img src="data:image/png;base64, <%= value %>" style="width:50px;height:50px;">
    </div>
    """

    # initiallizing the bokeh figure using the previous template for each functional group
    columns = [
        TableColumn(field="groups", title="Functional Groups"),
        TableColumn(field="images", title="Images", width=200, formatter=HTMLTemplateFormatter(template=template))
    ]
    num_groups = len(contained_functional_groups)

    table_height = min(200 + num_groups * 60, 800)

    data_table = DataTable(source=source, columns=columns, width=250, height=table_height, row_height=60)

    return data_table






