# -         MASSIVEChem       - 
 - Python package for applied analytical chemistry focused primarily on mass speectrometry 
#### Project in practical programming in chemistry course -- EPFL CH-200

## Package description 
MASSIVEChem, which stands for "Mass Analytical Spectrometry System for Investigation and Visual Extrapolation in Chemistry", is a pip-installable package developped at EPFL in 2024 focused on, as its name would suggest, analytical chemistry.
The aim of this package is to provide the user functions in order to simulate the mass spectrum of a molecule and to display this spectrum on a graph. The package also provides other features that can facilitate the chemical analysis of a molecule such as a functional group finder and an instauration calculator.

Developpers:
- Thomas Viking Christiansson, student in chemical engineering at EPFL    https://github.com/ThomasCsson
- Igor Gonteri, student in chemistry at EPFL                             https://github.com/igorgonteri
- Arthur Humery, student in chemical engineering at EPFL                https://github.com/Arthurhmy

### What is mass spectrometry ?
   - Mass spectrometry is an analytical technique used to identify and quantify chemical compounds in a sample by measuring the mass and sometimes the charge of molecules. It involves separating pre-charged ions according to their mass-to-charge ratio (m/z), then detecting and analysing them. This method is widely used in chemistry, biochemistry, pharmacology and other fields to characterise substances and understand their composition.

Now, let us go through the steps required to use this package !

## Installation

MASSIVEChem can be installed using pip as
```bash
pip install MASSIVEChem
```
The package can also be installed from source by running the following commands

First, clone the repository from github

```bash
git clone https://github.com/ThomasCsson/ppchem-project-Christiansson-Gonteri-Humery.git
cd ThomasCsson/ppchem-project-Christiansson-Gonteri-Humery
```
Then, execute the shell script. The shell scripts require two arguments, python version and gpu/cpu.

```bash
source scripts/install_deepchem_conda.sh 3.10 cpu
```

## Requirments
The package runs on python 3.10 but supports python 3.8 through 3.10
The package requires several other packages to function correctly.

```bash
matplotlib
bokeh
rdkit
pandas
```
If everything runs in order during the installation, the preceding packages should install automatically.
But check that the following packages are correctly installed using 

```bash
pip show "name of the package"
```

If not, install them using the following commands, otherwise the package will not work. 

```bash
pip install matplotlib
pip install bokeh
pip install rdkit
pip install pandas
```

## Usage

'''Show the most important function and use'''

## Getting started

To begin to use the package the following jupyter notebook will give you information about all the package's functions:

'''link to jupter notebook'''

