
from bokeh.layouts import column
from bokeh.models import Button, CustomJS, Div
from bokeh.plotting import show


def button_display(mol_weight, unsat, nbr_atoms, atom_list):
    #---------------------------------------------------------------------------------------------#
    '''
    button_display(mol_weight, unsat, nbr_atoms, atom_list)
    
    Input: molecular weight (int), unsaturation (int), number of atoms (int) and contained atoms (list) of the molecule.
    
    Output: bokeh layout of 4 buttons that displays the information

    '''
    #---------------------------------------------------------------------------------------------#

    # Change the list of atom to a string to display
    formatted_elements = [f"{part.split(':')[0].strip()} : {part.split(':')[1].strip()}" for part in atom_list]

    # Join all formatted elements into a single string with a space between them
    atom_string = ' , '.join(formatted_elements)

    # Create a Div element to display the information with adaptive size
    info_div = Div(text="", styles={'background-color': '#f0f0f0', 'padding': '10px', 'border': '2px solid #ccc', 'border-radius': '5px'})
    
    # Create Buttons
    button1 = Button(label="Molecular Mass", button_type="success")
    button2 = Button(label="Unsaturation", button_type="success")
    button3 = Button(label="Number of Atoms", button_type="success")
    button4 = Button(label="Types of Atoms", button_type="success")
    
    # JavaScript code to update the Div with the input information
    callback1 = CustomJS(args=dict(div=info_div, info=f'{mol_weight} g/mol'), code="""
        div.text = info;
        div.change.emit();
    """)
    if unsat == 1:
        info = '1 unsaturation'
    else:
        info = (f'{unsat} unsaturations')
    callback2 = CustomJS(args=dict(div=info_div, info = info), code="""
        div.text = info;
        div.change.emit();
    """)
    
    callback3 = CustomJS(args=dict(div=info_div, info=f'{nbr_atoms} atoms'), code="""
        div.text = info;
        div.change.emit();
    """)
    
    callback4 = CustomJS(args=dict(div=info_div, info=atom_string), code="""
        div.text = info;
        div.change.emit();
    """)
    
    # Attach the callbacks to the buttons' click events
    button1.js_on_click(callback1)
    button2.js_on_click(callback2)
    button3.js_on_click(callback3)
    button4.js_on_click(callback4)
    
    # Create the layout
    layout = column(button1, button2, button3, button4, info_div)
    
    return layout

show(button_display(230, 1, 13,['C : 10', 'H : 24']))


