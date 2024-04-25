import pandas as pd

# Read the .txt file into a DataFrame
df = pd.read_csv('/Users/thomaschristiansson/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Thomas/abundance.txt', sep='\t', header=None, names=['Atom', 'Mass', 'Percentage'])






# Extract the 'Mass' column as a list
atom_masses = df['Mass'].tolist()

# Print the list of atom masses
print(atom_masses)
# Print the DataFrame
print(df)