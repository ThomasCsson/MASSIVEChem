from bokeh.io import curdoc
from bokeh.layouts import column
from bokeh.models import TextInput, Button, Paragraph
from bokeh.plotting import figure

# Function to update the result
def update():
    # Get the input value and multiply by 2
    try:
        result = float(input_field.value) * 2
        result_paragraph.text = f"Result: {result}"
    except ValueError:
        result_paragraph.text = "Please enter a valid number."

# Create TextInput for user input
input_field = TextInput(value="1", title="Enter a number:")

# Create Button to trigger update
update_button = Button(label="Update")

# Bind update function to button click event
update_button.on_click(update)

# Create Paragraph to display result
result_paragraph = Paragraph(text="Result: ")

# Create a layout for the application
layout = column(input_field, update_button, result_paragraph)

# Add layout to current document
curdoc().add_root(layout)