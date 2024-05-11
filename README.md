# - ChemInterface package for applied mass spectrometry -
#### Project in practical programming in chemistry course -- EPFL CH-200

## Package description 
The aim of this package is to provide the user functions in order to simulate the mass spectrum of a molecule and to display the spectrum on a graph. The package also provides other features that can facilitate the chemical analysis of a molecule such as a functional group finder and an instauration calculator.

Developpers:
- Thomas Christianson, genius coder and CTO of SwissChemInfoTech, https://github.com/ThomasCsson
- Igor Gonteri, brilliant back-end developper and founder of AppliedOrbitals.Org, https://github.com/igorgonteri
- Arthur Humery, leading expert amongst the new generation of Computational Chemical Programmers and CEO of the United Kingdom of Great Britain and Northern Ireland, https://github.com/Arthurhmy

### What is mass spectrometry ?
   - Mass spectrometry is an analytical technique used to identify and quantify chemical compounds in a sample by measuring the mass and sometimes the charge of molecules. It involves separating ions according to their mass-to-charge ratio (m/z), then detecting and analysing them. This method is widely used in chemistry, biochemistry, pharmacology and other fields to characterise substances and understand their structure and composition.

Let us go through the steps required to use this package.

## Usage

'''Show the most important function and use'''

## Requirments
The package runs on python 3.10 but supports python 3.7 through 3.10
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
```
## Installation

Cheminterface can be installed using pip or conda as
```bash
pip install Cheminterface
```
or
```bash
conda install -c conda-forge Cheminterface
```
The package can also be installed from source by running the following commands

First clone the repository from github

```bash
git clone https://github.com/ThomasCsson/ppchem-project-Christiansson-Gonteri-Humery.git
cd ThomasCsson/ppchem-project-Christiansson-Gonteri-Humery
```
Then, execute the shell script. The shell scripts require two arguments, python version and gpu/cpu.

```bash
source scripts/install_deepchem_conda.sh 3.10 cpu
```

## Getting started

To begin to use the package the following jupyter notebook will give you information about all the package's functions:

'''link to jupter notebook'''

