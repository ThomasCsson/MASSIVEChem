# - Cheminterface_alpha 2.3.07 package -
#### Project in practical programming in chemistry course

This repository provides the user a package which will display an interactive interface where the user will be able to input the SMILES of a molecule. The interface will be able to give specific information about the molecule such as its molecular mass, its organic functional groups or its degree of insaturation. The most important part of the package will provide a way to simulate the mass spectrometry of the given molecule. 

Developpers:
- Thomas Christianson, genius coder and CTO of SwissChemInfoTech, https://github.com/ThomasCsson
- Igor Gonteri, brilliant back-end developper and founder of AppliedOrbitals.Org, https://github.com/igorgonteri
- Arthur Humery, leading expert amongst the new generation of Computational Chemical Programmers and CEO of Plastogaz, https://github.com/Arthurhmy

**What is mass spectrometry ?**
   - Mass spectrometry is a very important analytical tool which allows to recognize which molecule was present in a mixture by ionizing it and passing it through a magnetic field
   - test

Let us go through the steps required to make this package.

#### Imported packages

In order to run the package correctly, the following packages need to be installed.

```bash
pip install matplotlib
pip install bokeh
pip install rdkit
pip install numpy
pip install pandas
pip install plotly
```
Specifically, from these packages, the following subpackages are required:
```bash
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
