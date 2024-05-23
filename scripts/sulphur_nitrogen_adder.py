def sulphur_nitrogen_adder(x_in, y_in, has_N, has_S, count_N, count_S) -> list[float]:
    #---------------------------------------------------------------------------------------------#
    '''
    sulphur_nitrogen_adder(x_in, y_in, has_N, has_S, count_N, count_S)
    
    Input: two lists + 4 booleans:
    1. ordered list of the masses with values on y merged together if peaks within precision of apparatus
    2. ordered list of the probabilities of apparation of each of the molecules
    3. boolean that tells if the molecule contains Nitrogen
    4. boolean that tells if the molecule contains Sulphur
    5. number of Nitrogen atoms in the molecule
    6. number of Sulphur atoms in the molecule
    
    Output: two lists:
    1. ordered list of the masses with values on y merged together if peaks within precision of apparatus
    2. ordered list of the probabilities of apparation of each of the molecules
    
    (the mass in list 1 at index i is associated to the probability at index i in list 2)
    '''
    #---------------------------------------------------------------------------------------------#

    if not x_in:
        raise ValueError("Enter a non-empty list")
    if not y_in:
        raise ValueError("Enter a non-empty list")
    if len(x_in) != len(y_in):
        raise ValueError("The two lists must have the same length")
    if has_N not in [True, False]:
        raise ValueError("The third argument must be a boolean")
    if has_S not in [True, False]:
        raise ValueError("The fourth argument must be a boolean")
    if count_N < 0:
        raise ValueError("The fifth argument must be a positive integer")
    if count_S < 0:
        raise ValueError("The sixth argument must be a positive integer")
    if type(count_N) != int:
        raise ValueError("The fifth argument must be an integer")
    if type(count_S) != int:
        raise ValueError("The sixth argument must be an integer")
    if x_in != sorted(x_in):
        raise ValueError("The first list must be ordered")

    maximum = max(y_in)

    if has_N:
        x_in.append(x_in[1] - 0.006)  
        y_in.append(0.0035*count_N*maximum)  
    
    if has_S:
        x_in.append(x_in[1]-0.004)  
        y_in.append(0.008*count_S*maximum)


    return x_in, y_in
