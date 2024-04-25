import pandas as pd

df = pd.read_csv('/Users/thomaschristiansson/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Thomas/abundance.txt', sep='\t', header=None, names=['Atom', 'Mass', 'Percentage'])







mass = df['Mass'].tolist()
mass = [float(m) for m in mass]

abundance_percent = df['Percentage'].tolist()
abundance_percent = [float(ap) for ap in abundance_percent]

isotopes = df['Atom'].tolist()


print(mass)