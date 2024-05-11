from bokeh.io import show
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import DataTable, TableColumn
from bokeh.models import HTMLTemplateFormatter



def functional_group_display(groups_list):
    #---------------------------------------------------------------------------------------------#
    '''
    functional_group_display(groups_list)
    
    Input: list of the groups present in the molecule
    
    Output: Bokeh table with the names of the present functional groups as well as an image of each present functional group
    '''
    #---------------------------------------------------------------------------------------------#

    # dictionnary of the images of all functional groups

    functional_groups_images = {
        'Alcohol': 'Functional groups images/Alcohol_image.png',
        'Aldehyde': 'Functional groups images/Aldehyde_image.png',
        'Ketone': 'Functional groups images/Ketone_image.png',
        'Carboxylic Acid': 'Functional groups images/Acid_image.png',
        'Ester': 'Functional groups images/Ester_image.png',
        'Ether': 'Functional groups images/Ether_image.png',
        'Amide': 'Functional groups images/Amide_image.png',
        'Amine': 'Functional groups images/Amine_image.png',
        'Nitrile': 'Functional groups images/Nitrile_image.png',
        'Chloride': 'Functional groups images/Halogen_image.png',
        'Alkene': 'Functional groups images/Alkene_image.png',
        'Alkyne': 'Functional groups images/Alkyne_image.png',
        'Imine': 'Functional groups images/Imine_image.png',
        'Amino acid': 'Functional groups images/Amino_acid_image.png',
        'Proline': 'Functional groups images/Proline_image.png',
        'Thiol': 'Functional groups images/Thiol_image.png',
        'Sulfides': 'Functional groups images/Sulfides_image.png',
        'Acyl Chloride': 'Functional groups images/Acyl_chloride_image.png',
        'Anhydride': 'Functional groups images/Anhydride_image.png',
        'Nitro': 'Functional groups images/Nitro_image.png',
        'Enamine': 'Functional groups images/Enamine_image.png',
        'Enamine2': 'Functional groups images/Enamine2_image.png',
        'Enamine3': 'Functional groups images/Enamine3_image.png',
        'Imide': 'Functional groups images/Imide_image.png',
        'Azide': 'Functional groups images/Azide_image.png',
        'Enol': 'Functional groups images/Enol_image.png',
        'Hemiacetal': 'Functional groups images/Hemiacetal_image.png',
        'Carbonate': 'Functional groups images/Carbonate_image.png',
        'Carbonate2': 'Functional groups images/Carbonate2_image.png',
        'Disulfide': 'Functional groups images/Disulfide_image.png',
        'Sulfoxide': 'Functional groups images/Sulfoxide_image.png',
        'Sulfone': 'Functional groups images/Sulfone_image.png',
        'Sulfonic acid': 'Functional groups images/Sulfonic_acid_image.png',
        'Thioester': 'Functional groups images/Thioester_image.png',
        'Phosphine': 'Functional groups images/Phosphine_image.png',
        'Phosphate ester': 'Functional groups images/Phosphate_image.png',
        'Benzene': 'Functional groups images/Benzene_image.png',
        'Peroxide': 'Functional groups images/Peroxide_image.png'
}
    
    # creates a dictionnary of the present groups and associated images of the molecule

    present_group_images = []
    for x in groups_list:
        present_group_images.append(functional_groups_images[x])
    data = dict(
        groups=groups_list,
        images=[f'<img src="{group_image}" style="width:50px;height:50px;">' for group_image in present_group_images]
    )
    source = ColumnDataSource(data)

    #template for the bokeh table

    template = """
    <div>
    <%= value %>
    </div>
    """

    # initiallizing the bokeh figure using the previous template for each functional group

    columns = [
        TableColumn(field="groups", title="Functional Groups"),
        TableColumn(field="images", title="Images", width=200, formatter=HTMLTemplateFormatter(template=template))
    ]
    num_groups = len(groups_list)

    table_height = min(200 + num_groups * 60, 800)

    data_table = DataTable(source=source, columns=columns, width=250, height=table_height, row_height=60)

    return data_table



groups = ['Proline', 'Benzene', 'Aldehyde','Carboxylic Acid','Phosphine']  # Example with more than 3 groups


show(functional_group_display(groups))


