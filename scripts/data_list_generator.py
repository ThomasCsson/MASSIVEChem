import pandas as pd

def data_list_generator():

    #---------------------------------------------------------------------------------------------#
    '''
    data_list_generator ()

    Input: none

    Output: three lists:
    1. mass: [mass1, mass2, mass3,...]  
    2. abundance: [ab1, ab2, ab3,...]
    3. isotopes: [iso1, iso2, iso3,...]
    '''
    #---------------------------------------------------------------------------------------------#


    #Turn data of (Symbol | Mass | Probability) into lists 

    df = pd.read_csv('data/abundance.txt'
                    , sep='\t'
                    , header=None
                    , names=['Atom', 'Mass', 'Percentage'])


    #mass = [mass1, mass2, mass3,...]

    mass = df['Mass'].tolist()


    #change from elements to floats

    mass = [float(m) for m in mass]


    #abundance = [ab1, ab2, ab3,...]

    abundance_percent = df['Percentage'].tolist()


    #change from elements to floats

    abundance_percent = [float(ap) for ap in abundance_percent]


    #from percent to proba

    abundance = []
    for percent in abundance_percent:
        abundance.append(percent/100)


    #isotopes = [iso1, iso2, iso3,...]

    isotopes = df['Atom'].tolist()

    return mass, abundance, isotopes

print(data_list_generator())