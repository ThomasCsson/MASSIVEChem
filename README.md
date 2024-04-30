# - ChemInterface package for applied mass spectrometry -
#### Project in practical programming in chemistry course
#### EPFL CH-200

This repository provides the user a package which will display an interactive interface where the user will be able to input the SMILES of a molecule. The interface will be able to give specific information about the molecule such as its molecular mass, its organic functional groups or its degree of insaturation. The most important part of the package will provide a way to simulate the mass spectrometry of the given molecule. 

Developpers:
- Thomas Christianson, genius coder and CTO of SwissChemInfoTech, https://github.com/ThomasCsson
- Igor Gonteri, brilliant back-end developper and founder of AppliedOrbitals.Org, https://github.com/igorgonteri
- Arthur Humery, leading expert amongst the new generation of Computational Chemical Programmers and CEO of the United Kingdom of Great Britain and Northern Ireland, https://github.com/Arthurhmy

**What is mass spectrometry ?**
   - Mass spectrometry is an analytical technique used to identify and quantify chemical compounds in a sample by measuring the mass and sometimes the charge of molecules. It involves separating ions according to their mass-to-charge ratio (m/z), then detecting and analysing them. This method is widely used in chemistry, biochemistry, pharmacology and other fields to characterise substances and understand their structure and composition.

Let us go through the steps required to make this package.

#### Imported packages

In order to run the package correctly, the following packages need to be installed before running any pip-installed functions using the following command

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
