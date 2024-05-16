from bokeh.models import ColumnDataSource, HTMLTemplateFormatter
from bokeh.models.widgets import DataTable, TableColumn
from bokeh.plotting import show

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
        'Alcohol': '../data/Functional groups images/Alcohol_image.png',
        'Aldehyde': '../data/Functional groups images/Aldehyde_image.png',
        'Ketone': '../data/Functional groups images/Ketone_image.png',
        'Carboxylic Acid': '../data/Functional groups images/Acid_image.png',
        'Ester': '../data/Functional groups images/Ester_image.png',
        'Ether': '../data/Functional groups images/Ether_image.png',
        'Amide': '../data/Functional groups images/Amide_image.png',
        'Amine': '../data/Functional groups images/Amine_image.png',
        'Nitrile': '../data/Functional groups images/Nitrile_image.png',
        'Chloride': '../data/Functional groups images/Halogen_image.png',
        'Bromide': '../data/Functional groups images/Bromide_image.png',
        'Fluoride': '../data/Functional groups images/Fluoride_image.png',
        'Iodide': '../data/Functional groups images/Iodide_image.png',
        'Alkene': '../data/Functional groups images/Alkene_image.png',
        'Alkyne': '../data/Functional groups images/Alkyne_image.png',
        'Imine': '../data/Functional groups images/Imine_image.png',
        'Amino acid': '../data/Functional groups images/Amino_acid_image.png',
        'Proline': '../data/Functional groups images/Proline_image.png',
        'Thiol': '../data/Functional groups images/Thiol_image.png',
        'Sulfides': '../data/Functional groups images/Sulfides_image.png',
        'Acyl Chloride': '../data/Functional groups images/Acyl_chloride_image.png',
        'Anhydride': '../data/Functional groups images/Anhydride_image.png',
        'Nitro': '../data/Functional groups images/Nitro_image.png',
        'Enamine': '../data/Functional groups images/Enamine_image.png',
        'Enamine2': '../data/Functional groups images/Enamine2_image.png',
        'Enamine3': '../data/Functional groups images/Enamine3_image.png',
        'Imide': '../data/Functional groups images/Imide_image.png',
        'Azide': '../data/Functional groups images/Azide_image.png',
        'Enol': '../data/Functional groups images/Enol_image.png',
        'Hemiacetal': '../data/Functional groups images/Hemiacetal_image.png',
        'Carbonate': '../data/Functional groups images/Carbonate_image.png',
        'Carbonate2': '../data/Functional groups images/Carbonate2_image.png',
        'Disulfide': '../data/Functional groups images/Disulfide_image.png',
        'Sulfoxide': '../data/Functional groups images/Sulfoxide_image.png',
        'Sulfone': '../data/Functional groups images/Sulfone_image.png',
        'Sulfonic acid': '../data/Functional groups images/Sulfonic_acid_image.png',
        'Thioester': '../data/Functional groups images/Thioester_image.png',
        'Phosphine': '../data/Functional groups images/Phosphine_image.png',
        'Phosphate ester': '../data/Functional groups images/Phosphate_image.png',
        'Benzene': '../data/Functional groups images/Benzene_image.png',
        'Peroxide': '../data/Functional groups images/Peroxide_image.png'
}
    
    # creates a dictionnary of the present groups and associated images of the molecule

    present_group_images = []

    for x in groups_list:
        if x in functional_groups_images.keys():
            present_group_images.append(functional_groups_images[x])
        else:
            pass

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

show(functional_group_display(['Alcohol', 'InvalidGroup']))

import unittest

class TestFunctionalGroupDisplay(unittest.TestCase):
    def test_functional_group_display(self):
        # Test the function with a list of functional groups
        groups_list = ['Alcohol', 'Aldehyde', 'Ketone']
        table = functional_group_display(groups_list)
        self.assertIsNotNone(table)
        self.assertEqual(len(table.columns), 2)  # Check if there are two columns
        self.assertEqual(table.columns[0].field, "groups")  # Check if the first column is for functional groups
        self.assertEqual(table.columns[1].field, "images")  # Check if the second column is for images

    def test_functional_group_display_empty_list(self):
        # Test the function with an empty list
        table = functional_group_display([])
        self.assertIsNotNone(table)
        self.assertEqual(len(table.columns), 2)  # Check if there are two columns
        self.assertEqual(len(table.source.data['groups']), 0)  # Check if the number of groups is 0

    def test_functional_group_display_invalid_group(self):
        # Test the function with an invalid group
        groups_list = ['Alcohol', 'InvalidGroup']
        table = functional_group_display(groups_list)
        self.assertIsNotNone(table)
        self.assertEqual(len(table.columns), 2)  # Check if there are two columns
        self.assertEqual(len(table.source.data['groups']), 1)  # Check if the number of groups is 1
        self.assertEqual(table.source.data['groups'][0], 'Alcohol')  # Check if the valid group is displayed

if __name__ == '__main__':
    unittest.main()

