import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from mordred import Calculator, descriptors
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def calculate_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    calc = Calculator(descriptors, ignore_3D=True)
    descriptors = calc(mol)
    descriptor_values = [d.value for d in descriptors]
    return descriptor_values

def predict_nmr_spectrum(smiles):
    descriptor_values = calculate_descriptors(smiles)
    # Here, you would typically use a machine learning model to predict the NMR spectrum
    # For demonstration purposes, let's just generate a random spectrum
    nmr_spectrum = np.random.rand(100)
    return nmr_spectrum

def plot_nmr_spectrum(nmr_spectrum):
    plt.figure(figsize=(10, 6))
    plt.plot(nmr_spectrum)
    plt.title("Predicted NMR Spectrum")
    plt.xlabel("Chemical Shift")
    plt.ylabel("Intensity")
    plt.grid(True)
    plt.show()

# Example SMILES string
smiles = "CCO"
# Predict NMR spectrum
nmr_spectrum = predict_nmr_spectrum(smiles)
# Plot NMR spectrum
plot_nmr_spectrum(nmr_spectrum)

