![logo](IMG_2856.jpg)
[![GitHub](https://img.shields.io/badge/github-%23121011.svg?style=for-the-badge&logo=github&logoColor=white)](https://github.com/ThomasCsson/MASSIVEChem)
[![Thomas3](https://img.shields.io/badge/Python-FFD43B?style=for-the-badge&logo=python&logoColor=blue)](https://www.python.org/)
![Thomas4](https://img.shields.io/badge/HTML5-E34F26?style=for-the-badge&logo=html5&logoColor=white)
[![Thomas6](https://img.shields.io/badge/Jupyter-F37626.svg?&style=for-the-badge&logo=Jupyter&logoColor=purple)](https://jupyter.org/)


# -         MASSIVEChem       - 

[![GitHub2](https://img.shields.io/badge/Maintained%3F-yes-turquoise.svg)](https://pypi.org/user/Arthur.hmy/)
[![python](https://img.shields.io/badge/Python-3.10-3776AB.svg?style=flat&logo=python&logoColor=orange)](https://www.python.org)
[![GitHub3](https://img.shields.io/badge/Contributors-3-green.svg)](https://github.com/ThomasCsson/MASSIVEChem/graphs/contributors)
[![GitHub3](https://img.shields.io/badge/License-3-purple.svg)](https://github.com/ThomasCsson/MASSIVEChem/blob/main/LICENSE.txt)
[![GitHub3](https://img.shields.io/badge/EPFL-CH200-red.svg)](https://edu.epfl.ch/studyplan/en/bachelor/chemistry-and-chemical-engineering/coursebook/practical-programming-in-chemistry-CH-200)


 - Python package for applied analytical chemistry focused primarily on mass speectrometry 
#### Project within _practical programming in chemistry_ course -- EPFL CH-200

## Package description
[![GitHub3](http://ForTheBadge.com/images/badges/built-with-science.svg)](https://x.com/pschwllr/status/1760713822111723990)
[![GitHub4](http://ForTheBadge.com/images/badges/made-with-python.svg)](https://www.python.org/)

MASSIVEChem, which stands for "Mass Analytical Spectrometry System for Investigation and Visual Extrapolation in Chemistry", is a pip-installable python package developped at EPFL in 2024 focused on, as its name would suggest, analytical chemistry.
The aim of this package is to provide the user with functions in order to simulate the mass spectrum of a molecule and to display this spectrum on a graph. The package also provides other features that can facilitate the chemical analysis of a molecule such as a functional group finder and an unsaturation calculator.

Developpers:
- Thomas Viking Christiansson, student in chemical engineering at EPFL    [![jhc github](https://img.shields.io/badge/GitHub-ThomasCsson-181717.svg?style=flat&logo=github)](https://github.com/ThomasCsson)
- Igor Gonteri, student in chemistry at EPFL                             [![jhc github](https://img.shields.io/badge/GitHub-igorgonteri-181717.svg?style=flat&logo=github)](https://github.com/igorgonteri)
- Arthur Humery, student in chemical engineering at EPFL                [![jhc github](https://img.shields.io/badge/GitHub-Arthurhmy-181717.svg?style=flat&logo=github)](https://github.com/Arthurhmy)

### What is mass spectrometry ?
Mass spectrometry (MS) is an analytical technique used to measure the mass-to-charge ratio of ions. It helps identify the structure of the chemical compound present in a sample by generating a spectrum of the masses of its ions. The process involves three main steps:

Ionization: The sample is ionized, which means its molecules are converted into charged particles (ions). This can be done using various methods, such as electron impact (EI), electrospray ionization (ESI), or matrix-assisted laser desorption/ionization (MALDI). In this case, the ionisation is set to the most commonly used method; deprotonation.

Mass Analysis: The ions are separated based on their mass-to-charge ratio (m/z). This is usually done using a mass analyzer, such as a quadrupole, time-of-flight (TOF), or an ion trap. Each type of mass analyzer works differently but ultimately serves to distinguish ions by their specific m/z values. The unit of these m/z values is 1 Th or 1 $\frac{Da}{e}$.

On the y axis of the output spectrum, the relative abundance of the different ions is plotted. This abundance is given by the different natural abundances of the different isotopes of the atoms in teh molecule. For example, the relative abundance of $^{13}C$ is of 1.1% and that of $^{12}C$ is of 98.9%.

Mass spectrometry is widely used in various fields, including:

Chemistry: For molecular identification and structural elucidation.
Biochemistry: For studying proteins, peptides, and other biomolecules.
Pharmaceuticals: For drug development and metabolite analysis.
Environmental Science: For detecting pollutants and analyzing environmental samples.
Clinical Diagnostics: For identifying biomarkers and analyzing complex biological samples.
The technique's sensitivity, accuracy, and ability to analyze complex mixtures make it an essential tool in both research and applied sciences.

Now, let us go through the steps required to use this package!

## Installation
[![Thomas2](https://img.shields.io/badge/pypi-3775A9?style=for-the-badge&logo=pypi&logoColor=white)](https://pypi.org/project/MASSIVEChem/)

MASSIVEChem can be installed using pip as
```bash
pip install MASSIVEChem
```


[![GitHub](https://img.shields.io/badge/github-%23121011.svg?style=for-the-badge&logo=github&logoColor=white)](https://github.com/ThomasCsson/MASSIVEChem?tab=readme-ov-file)

Alternatively, the package can be directly installed from the GitHub repository via pip using the following command in the terminal
```bash
pip install git+https://github.com/ThomasCsson/MASSIVEChem
```


![Thomas1](https://img.shields.io/badge/GIT-E44C30?style=for-the-badge&logo=git&logoColor=white)

The package can also be installed from source by running the following commands

First, clone the repository from github and go in the folder. 
```bash
git clone https://github.com/ThomasCsson/MASSIVEChem.git
cd MASSIVEChem
```
Then, install the package using : 
```bash
pip install -e .
```







## Requirments
The package runs on python 3.10 but supports python 3.8 through 3.10.
It requires several other packages to function correctly.

```bash
matplotlib
bokeh
rdkit
pandas
panel
```

(Note: panel is not directly used, but is required for added functionality involving the 3 dimmensinal visualisation of the input molecule)

If all goes well during installation, the preceding packages should all install automatically.
But this can be checked by veryfying that they have all been installed in the desired environment. To do this, simply write the following command in the terminal:

```bash
pip show "name of the package"
```

If not, install them using the following commands. (Bear in mind that the package will not run without its dependencies. 

```bash
pip install matplotlib
pip install bokeh
pip install rdkit
pip install pandas
pip install panel
```
Additionally, the package 'xyz2graph' is required to run the 3D imaging functionallity in the function XXX. This package is not pip-installable, so to intall it yourslef, the following command needs to be run.

```bash
python -m pip install git+https://github.com/zotko/xyz2graph.git
```


## Usage

The principal function of this package takes the SMILEs of a molecule as an input and displays the mass spectrometry of the molecule as well as the molecule itself and  the functional groups it contains.

An example on how to make the function work is shown below for benzylpenicilin:

The ionization method is set to monodeprotonation and the resolution of the apparatus is 0.01 Th

```bash
import MASSIVEChem.MASSIVEChem as ms
from bokeh.plotting import show

mol_smi = 'CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C'
apparatus_resolution = 0.01

show(ms.spectrum(mol_smi, True, apparatus_resolution))

'''
The first input in ms.spectrum is the molecule under SMILEs representation,
the second computes an approximate spectrum if True and the precise spectrum if False
and the third is the resolution of the apparatus (typically, this value is of 0.01
```

The output of this command will be:

![Spectrum](Spectrum_output.png)

Note that here there appear to be two overlapping peaks at ~ 334 [th]. This is due to the presence of an odd number of sulphur atoms. This causes there to be a peak 0.004 [th] infront of the second peak. This can be verified by zooming in on the sectrum:

<img width="600" alt="Focused_spectrum" src="https://github.com/ThomasCsson/MASSIVEChem/assets/160872481/440ac2ea-c2fe-40ff-b1ed-1a447024bcb0">

## Getting started
[![jupyter](https://img.shields.io/badge/Jupyter-Lab-F37626.svg?style=flat&logo=Jupyter)](https://jupyterlab.readthedocs.io/en/stable)

To begin to use the package the following jupyter notebook will give you information about all the package's functions:

'''link to jupter notebook'''

