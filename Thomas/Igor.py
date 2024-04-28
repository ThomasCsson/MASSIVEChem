
#importing packages
import matplotlib.pyplot as plt
import pandas as pd


# Lecture du fichier texte
with open('abundance.txt', 'r', encoding='ISO-8859-1') as file:
    # Sauter la première ligne
    next(file)
    lines = file.readlines()

# Initialisation des listes pour stocker les valeurs des colonnes
colonne1 = []
colonne2 = []
colonne3 = []


# Parcours des lignes et extraction des valeurs
for line in lines:

    values = line.split()

    colonne1.append(values[0])
    colonne2.append(float(values[1]))
    colonne3.append(float(values[2]))

data = {
    'Atom': colonne1,
    'Mass': colonne2,
    'Abundance %': pd.to_numeric(colonne3)
}

# Création du DataFrame
DF = pd.DataFrame(data)

print(DF)







import pandas as pd

# Read the .txt file into a DataFrame
df = pd.read_csv('atom_masses.txt', sep='\t', header=None, names=['Atom', 'Mass', 'Percentage'])

# Print the DataFrame
print(df)