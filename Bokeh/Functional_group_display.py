from bokeh.io import show
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import DataTable, TableColumn
from bokeh.models import HTMLTemplateFormatter

groups = ['Sulfone', 'Benzene', 'Imide', 'Alcohol', 'Aldehyde', 'Ketone']  # Example with more than 3 groups
functional_groups_images = {
    'Alcohol': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Alcohol_image.png',
    'Aldehyde': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Aldehyde_image.png',
    'Ketone': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Ketone_image.png',
    'Carboxylic Acid': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Acid_image.png',
    'Ester': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Ester_image.png',
    'Ether': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Ether_image.png',
    'Amide': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Amide_image.png',
    'Amine': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Amine_image.png',
    'Nitrile': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Nitrile_image.png',
    'Chloride': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Halogen_image.png',
    'Alkene': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Alkene_image.png',
    'Alkyne': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Alkyne_image.png',
    'Imine': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Imine_image.png',
    'Amino acid': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Amino_acid_image.png',
    'Thiol': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Thiol_image.png',
    'Sulfides': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Sulfides_image.png',
    'Acyl Chloride': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Acyl_chloride_image.png',
    'Anhydride': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Anhydride_image.png',
    'Nitro': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Nitro_image.png',
    'Enamine': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Enamine_image.png',
    'Imide': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Imide_image.png',
    'Azide': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Azide_image.png',
    'Enol': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Enol_image.png',
    'Hemiacetal': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Hemiacetal_image.png',
    'Carbonate': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Carbonate_image.png',
    'Disulfide': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Disulfide_image.png',
    'Sulfoxide': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Sulfoxide_image.png',
    'Sulfone': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Sulfone_image.png',
    'Sulfonic acid': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Sulfonic_acid_image.png',
    'Thioester': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Thioester_image.png',
    'Phosphine': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Phosphine_image.png',
    'Phosphate ester': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Phosphate_image.png',
    'Benzene': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Benzene_image.png'
}

def functional_group_display(groups):
    functional_groups_images = {
        'Alcohol': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Alcohol_image.png',
        'Aldehyde': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Aldehyde_image.png',
        'Ketone': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Ketone_image.png',
        'Carboxylic Acid': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Acid_image.png',
        'Ester': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Ester_image.png',
        'Ether': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Ether_image.png',
        'Amide': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Amide_image.png',
        'Amine': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Amine_image.png',
        'Nitrile': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Nitrile_image.png',
        'Chloride': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Halogen_image.png',
        'Alkene': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Alkene_image.png',
        'Alkyne': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Alkyne_image.png',
        'Imine': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Imine_image.png',
        'Amino acid': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Amino_acid_image.png',
        'Thiol': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Thiol_image.png',
        'Sulfides': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Sulfides_image.png',
        'Acyl Chloride': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Acyl_chloride_image.png',
        'Anhydride': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Anhydride_image.png',
        'Nitro': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Nitro_image.png',
        'Enamine': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Enamine_image.png',
        'Imide': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Imide_image.png',
        'Azide': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Azide_image.png',
        'Enol': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Enol_image.png',
        'Hemiacetal': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Hemiacetal_image.png',
        'Carbonate': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Carbonate_image.png',
        'Disulfide': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Disulfide_image.png',
        'Sulfoxide': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Sulfoxide_image.png',
        'Sulfone': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Sulfone_image.png',
        'Sulfonic acid': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Sulfonic_acid_image.png',
        'Thioester': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Thioester_image.png',
        'Phosphine': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Phosphine_image.png',
        'Phosphate ester': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Phosphate_image.png',
        'Benzene': '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/Functional groups images/Benzene_image.png'
    }
    present_group_images = []
    for x in groups:
        present_group_images.append(functional_groups_images[x])
    data = dict(
        groups=groups,
        images=[f'<img src="{group_image}" style="width:50px;height:50px;">' for group_image in present_group_images]
    )
    source = ColumnDataSource(data)

    template = """
    <div>
    <%= value %>
    </div>
    """

    columns = [
        TableColumn(field="groups", title="Functional Groups"),
        TableColumn(field="images", title="Images", width=200, formatter=HTMLTemplateFormatter(template=template))
    ]
    num_groups = len(groups)

    table_height = min(200 + num_groups * 60, 800)

    data_table = DataTable(source=source, columns=columns, width=250, height=table_height, row_height=60)

    return data_table

show(functional_group_display(groups))


