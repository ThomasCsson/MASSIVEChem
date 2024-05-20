
import py3Dmol

# Create a Py3Dmol viewer
view = py3Dmol.view(width=400, height=400)

# Add a simple model to the viewer
view.addModel('water')

# Display the viewer
view.show()
