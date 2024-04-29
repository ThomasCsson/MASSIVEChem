import tkinter as tk

# Create a new Tkinter window
window = tk.Tk()

# Set the title of the window
window.title("Hello Tkinter")

# Create a label widget
label = tk.Label(window, text="Hello, Tkinter!")

# Pack the label widget into the window
label.pack()

# Start the Tkinter event loop
window.mainloop()
