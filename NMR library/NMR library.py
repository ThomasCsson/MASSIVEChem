import pandas as pd
import numpy as np

with open('Funct_groups_shift.txt', 'r', encoding='ISO-8859-1') as file:
    # Sauter la première ligne
    next(file)
    lines = file.readlines()

colonne1 = []
colonne2 = []
colonne3 = []
colonne4 = []

# Parcours des lignes et extraction des valeurs
for line in lines:
    values = line.split()  # Séparation des valeurs par espace

    # Vérifier si la ligne contient le nombre attendu de valeurs
    if len(values) >= 2:
        # Remplacement des virgules par des points pour les décimales et conversion en float
        colonne1.append(values[0])
        colonne2.append(float(values[1].replace(',', '.')))
        colonne3.append(float(values[2].replace(',', '.')))
        colonne4.append(float(values[3].replace(',', '.')))

Functional_group = {
    'Functional Group': colonne1,
    'Shift for M = Methyl': colonne2,
    'Shift for M = Methylene': colonne3,
    'Shift for M = Methile': colonne4
}

Functional_group_DF = pd.DataFrame(Functional_group)
print(Functional_group_DF)

def range1(x,y):
    result = [x + k * 0.05 for k in range(int((y - x) / 0.05) + 1)]
    return result

Ranges_functional_groups = {
    'Phenyl': range1(4, 10),
    'Alcohols': range1(1, 6),
    'Arenes': range1(9,6),
    'Alkenes': range1(4,8),
    'Ethers': range1(3,5),
    'Ketone': range1(2,3),
    'alkane': range1(1,2),
    'Methyl next to alkene': range1(1,2),
    'Methyl next to halogen': range1(1,2),

}
specific_func_group_shift = {
    'C1CC1': 0.5,
    'c1ccCcc1': 2.5,
    'CNC': 2.8,
    'CSC': 2.5,
    'CS': 3.8,
    'CN': 4.2,
    'C(=O)O': 11,
    'C=O': 10,
}