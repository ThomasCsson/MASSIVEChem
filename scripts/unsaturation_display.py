
from bokeh.layouts import column
from bokeh.models import Button, CustomJS, Div
from bokeh.plotting import show, output_notebook


def button_display(info1, info2, info3, info4):
    # Create a Div element to display the information with adaptive size
    info_div = Div(text="", styles={'background-color': '#f0f0f0', 'padding': '10px', 'border': '2px solid #ccc', 'border-radius': '5px'})
    
    # Create Buttons
    button1 = Button(label="Molecular Mass", button_type="success")
    button2 = Button(label="Unsaturation", button_type="success")
    button3 = Button(label="Number of Atoms", button_type="success")
    button4 = Button(label="Types of Atoms", button_type="success")
    
    # JavaScript code to update the Div with the input information
    callback1 = CustomJS(args=dict(div=info_div, info=f'{info1} g/mol'), code="""
        div.text = info;
        div.change.emit();
    """)
    if info2 == 1:
        info = '1 unsaturation'
    else:
        info = (f'{info2} unsaturations')
    callback2 = CustomJS(args=dict(div=info_div, info = info), code="""
        div.text = info;
        div.change.emit();
    """)
    
    callback3 = CustomJS(args=dict(div=info_div, info=f'{info3} atoms'), code="""
        div.text = info;
        div.change.emit();
    """)
    
    callback4 = CustomJS(args=dict(div=info_div, info=info4), code="""
        div.text = info;
        div.change.emit();
    """)
    
    # Attach the callbacks to the buttons' click events
    button1.js_on_click(callback1)
    button2.js_on_click(callback2)
    button3.js_on_click(callback3)
    button4.js_on_click(callback4)
    
    # Create the layout and display it
    layout = column(button1, button2, button3, button4, info_div)
    
    return layout

show(button_display(230, 1, 13, "C=3 N=2 O=2 H=6vvvvvvvvvvvr4nurtcvhineub fberikvniv bzinebz"))


