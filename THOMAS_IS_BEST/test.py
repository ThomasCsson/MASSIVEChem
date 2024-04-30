
list_atoms = ['C','C','N','N','N']

has_N = False
count_N = 0
has_S = False
count_S = 0
if 'N' in list_atoms and list_atoms.count('N')%2 == 1:
    has_N = True
    count_N = list_atoms.count('N')
elif 'S' in list_atoms and list_atoms.count('S')%2 == 1:
    has_S = True
    count_S = list_atoms.count('S')

print(has_N,has_S,count_N,count_S)