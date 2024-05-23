def peak_sorter(x_in, y_in) -> list[float]:
    #---------------------------------------------------------------------------------------------#
    '''
    peak_sorter(x_in, y_in)
    
    Input: two lists:
    1. ordered list of the masses with values on y merged together if peaks within precision of apparatus
    2. ordered list of the probabilities of apparation of each of the molecules
    
    Output: two lists:
    1. ordered list of the masses with values on y merged together if peaks within precision of apparatus
    2. ordered list of the probabilities of apparation of each of the molecules
    
    (the mass in list 1 at index i is associated to the probability at index i in list 2)
    '''
    #---------------------------------------------------------------------------------------------#

    x_out, y_out = [],[]
    while len(x_in)>0:
        min_x = min(x_in)
        index_min = x_in.index(min_x)
        x_out.append(min_x)
        y_out.append(y_in[index_min])
        x_in.pop(index_min)
        y_in.pop(index_min)

    return x_out, y_out