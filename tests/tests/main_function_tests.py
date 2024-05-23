def main_function(list_atoms, imprecision_True_False) -> list[float]:
    #---------------------------------------------------------------------------------------------#
    '''
    main_function(list_atoms)
    
    Input: list of atoms that enter the apparatus
    
    Output: two lists:
    1. list of the masses (of individual molecules) of each possible combination of isotopes
    2. list of the probabilities of apparation of each of the molecules 
    (the mass in list 1 at index i is associated to the probability at index i in list 2)
    '''
    #---------------------------------------------------------------------------------------------#
    render_imprecise_list = imprecision_True_False 
    #Set arg to be True for long molecules, set arg to False for short molecules/if precision for minuscule peaks is important
    
    mass = [108.904755, 106.90509, 26.981539, 39.962383, 37.96273, 35.967545, 74.92159, 196.96654, 11.009305, 10.012937, 137.90523, 136.9058, 135.90456, 134.90567, 133.90448, 131.90504, 129.90628, 9.012182, 208.98038, 78.918335, 80.91629, 13.003355, 12.0, 45.95369, 43.95548, 42.958767, 41.95862, 39.96259, 47.952534, 110.90418, 109.90301, 107.90418, 115.904755, 113.90336, 112.9044, 111.902756, 141.90924, 139.90543, 137.90599, 135.90714, 36.965904, 34.968853, 58.933197, 53.93888, 52.94065, 51.94051, 49.946045, 132.90543, 62.939598, 157.92441, 155.92528, 163.92917, 162.92873, 161.92679, 160.92693, 159.92519, 169.93546, 167.93237, 166.93205, 165.93028, 163.9292, 161.92877, 152.92122, 150.91985, 18.998404, 57.933277, 56.935394, 53.939613, 70.9247, 68.92558, 157.9241, 156.92395, 155.92212, 154.92262, 153.92087, 151.91978, 159.92705, 75.9214, 73.92118, 72.92346, 71.92208, 69.92425, 2.014, 1.007825, 4.0026, 3.01603, 173.94005, 179.94655, 178.94582, 177.9437, 176.94322, 203.97346, 201.97061, 200.97028, 199.9683, 198.96825, 197.96674, 195.9658, 164.93031, 126.90447, 114.90388, 112.90406, 190.96059, 192.96292, 40.961826, 39.964, 38.963707, 77.9204, 85.910614, 83.91151, 82.91414, 81.91348, 79.91638, 138.90634, 137.9071, 7.016003, 6.015121, 174.94077, 25.982594, 24.985838, 23.985043, 54.938046, 94.90584, 93.90508, 91.90681, 99.90748, 97.9054, 96.90602, 95.90468, 15.000108, 14.003074, 22.989767, 92.90638, 142.9098, 141.90771, 149.92088, 147.91689, 145.91312, 144.91257, 143.91008, 21.991383, 20.993843, 19.992435, 61.928345, 60.931057, 59.930786, 57.935345, 63.927967, 237.0482, 17.99916, 16.99913, 15.994915, 189.95844, 188.95813, 187.95586, 186.95573, 185.95383, 183.95248, 191.96147, 30.973763, 206.97588, 205.97444, 203.97302, 207.97662, 109.90517, 107.90389, 105.90348, 104.90508, 103.90403, 101.90563, 140.90765, 189.95992, 197.96786, 195.96492, 194.96477, 193.96266, 191.96101, 86.90919, 84.9118, 186.95575, 184.95296, 102.9055, 103.905426, 101.90435, 100.90558, 99.90422, 98.90594, 97.90529, 95.9076, 35.96708, 33.967865, 32.971455, 31.97207, 122.90421, 120.903824, 44.95591, 77.917305, 76.919914, 75.91921, 73.92248, 81.916695, 79.91652, 29.97377, 28.976496, 27.976927, 153.92221, 151.91972, 149.91727, 148.91718, 147.91483, 146.9149, 143.912, 123.90527, 121.90344, 119.9022, 118.90331, 117.90161, 116.902954, 115.90175, 114.90335, 113.90279, 111.90482, 87.90562, 86.90888, 85.90926, 83.91343, 180.948, 179.94746, 158.92534, 125.90331, 124.904434, 123.902824, 122.904274, 121.90305, 119.904045, 129.90623, 127.904465, 232.03806, 46.951763, 45.95263, 49.944794, 48.947872, 47.94795, 204.9744, 202.97232, 168.93422, 50.943962, 49.947163, 185.95436, 183.95093, 182.95023, 181.9482, 179.9467, 173.93886, 172.9382, 171.93637, 170.93633, 169.93475, 167.9339, 175.94257, 69.92532, 67.92484, 66.92713, 65.92603, 63.929146, 93.90647, 91.90504, 90.90565, 89.9047, 95.90827]
    abundance = [0.48161000000000004, 0.51839, 1.0, 0.996, 0.00063, 0.00337, 1.0, 1.0, 0.8009999999999999, 0.19899999999999998, 0.7170000000000001, 0.11230000000000001, 0.0785, 0.06593, 0.0242, 0.00101, 0.00106, 1.0, 1.0, 0.5069, 0.49310000000000004, 0.011000000000000001, 0.9890000000000001, 4e-05, 0.02086, 0.00135, 0.00647, 0.96941, 0.00187, 0.128, 0.1249, 0.0089, 0.07490000000000001, 0.2873, 0.1222, 0.2413, 0.11130000000000001, 0.8843000000000001, 0.0025, 0.0019, 0.24230000000000002, 0.7576999999999999, 1.0, 0.02365, 0.095, 0.8379000000000001, 0.043449999999999996, 1.0, 0.6917, 0.001, 0.0006, 0.282, 0.249, 0.255, 0.18899999999999997, 0.023399999999999997, 0.149, 0.268, 0.22949999999999998, 0.336, 0.0161, 0.0014000000000000002, 0.522, 0.478, 1.0, 0.0028000000000000004, 0.021, 0.059000000000000004, 0.39892000000000005, 0.60108, 0.2484, 0.1565, 0.2047, 0.14800000000000002, 0.0218, 0.002, 0.2186, 0.07440000000000001, 0.3594, 0.07719999999999999, 0.2766, 0.21239999999999998, 0.00015, 0.99985, 1.0, 1.37e-08, 0.0016200000000000001, 0.35100000000000003, 0.13629, 0.27297, 0.18606, 0.0687, 0.2986, 0.1318, 0.231, 0.16870000000000002, 0.09970000000000001, 0.0015, 1.0, 1.0, 0.9570000000000001, 0.043, 0.373, 0.627, 0.067302, 0.000117, 0.932581, 0.0034999999999999996, 0.17300000000000001, 0.57, 0.115, 0.11599999999999999, 0.0225, 0.999088, 0.000902, 0.925, 0.075, 0.9741, 0.1101, 0.1, 0.7898999999999999, 1.0, 0.1592, 0.0925, 0.1484, 0.09630000000000001, 0.2413, 0.0955, 0.1668, 0.0037, 0.9963, 1.0, 1.0, 0.12179999999999999, 0.2713, 0.0564, 0.0576, 0.17190000000000003, 0.083, 0.23800000000000002, 0.0925, 0.0027, 0.9048, 0.03634, 0.011399999999999999, 0.26222999999999996, 0.68077, 0.009260000000000001, 1.0, 0.002, 0.0004, 0.9976, 0.264, 0.161, 0.133, 0.016, 0.0158, 0.0002, 0.41, 1.0, 0.221, 0.24100000000000002, 0.013999999999999999, 0.524, 0.11720000000000001, 0.2646, 0.2733, 0.22329999999999997, 0.1114, 0.0102, 1.0, 0.0001, 0.07200000000000001, 0.253, 0.33799999999999997, 0.32899999999999996, 0.0079, 0.2783, 0.7217, 0.626, 0.374, 1.0, 0.18600000000000003, 0.316, 0.171, 0.126, 0.127, 0.018600000000000002, 0.0554, 0.0002, 0.0421, 0.0075, 0.9501999999999999, 0.4264, 0.5736, 1.0, 0.2377, 0.07629999999999999, 0.09359999999999999, 0.0089, 0.0874, 0.4961, 0.031, 0.0467, 0.9223, 0.22699999999999998, 0.267, 0.07400000000000001, 0.138, 0.113, 0.15, 0.031, 0.0579, 0.0463, 0.3259, 0.0858, 0.2422, 0.0768, 0.14529999999999998, 0.0036, 0.006500000000000001, 0.0097, 0.8258, 0.07, 0.0986, 0.005600000000000001, 0.9999800000000001, 0.00012, 1.0, 0.1893, 0.0712, 0.0479, 0.00905, 0.0259, 0.00095, 0.3387, 0.317, 1.0, 0.073, 0.08, 0.054000000000000006, 0.055, 0.738, 0.7047599999999999, 0.29524, 1.0, 0.9975, 0.0025, 0.28600000000000003, 0.307, 0.14279999999999998, 0.263, 0.0012, 0.318, 0.1612, 0.21899999999999997, 0.14300000000000002, 0.0305, 0.0013, 0.127, 0.006, 0.188, 0.040999999999999995, 0.27899999999999997, 0.486, 0.17379999999999998, 0.17149999999999999, 0.11220000000000001, 0.5145000000000001, 0.027999999999999997]
    isotopes = ['Ag', 'Ag', 'Al', 'Ar', 'Ar', 'Ar', 'As', 'Au', 'B', 'B', 'Ba', 'Ba', 'Ba', 'Ba', 'Ba', 'Ba', 'Ba', 'Be', 'Bi', 'Br', 'Br', 'C', 'C', 'Ca', 'Ca', 'Ca', 'Ca', 'Ca', 'Ca', 'Cd', 'Cd', 'Cd', 'Cd', 'Cd', 'Cd', 'Cd', 'Ce', 'Ce', 'Ce', 'Ce', 'Cl', 'Cl', 'Co', 'Cr', 'Cr', 'Cr', 'Cr', 'Cs', 'Cu', 'Dy', 'Dy', 'Dy', 'Dy', 'Dy', 'Dy', 'Dy', 'Er', 'Er', 'Er', 'Er', 'Er', 'Er', 'Eu', 'Eu', 'Fe ', 'Fe', 'Fe', 'Fe', 'Ga', 'Ga', 'Gd', 'Gd', 'Gd', 'Gd', 'Gd', 'Gd', 'Gd', 'Ge', 'Ge', 'Ge', 'Ge', 'Ge', 'H', 'H', 'He', 'He', 'Hf', 'Hf', 'Hf', 'Hf', 'Hf', 'Hg', 'Hg', 'Hg', 'Hg', 'Hg', 'Hg', 'Hg', 'Ho', 'I', 'In', 'In', 'Ir', 'Ir', 'K', 'K', 'K', 'Kr', 'Kr', 'Kr', 'Kr', 'Kr', 'Kr', 'La', 'La', 'Li', 'Li', 'Lu', 'Mg', 'Mg', 'Mg', 'Mn', 'Mo', 'Mo', 'Mo', 'Mo', 'Mo', 'Mo', 'Mo', 'N', 'N', 'Na', 'Nb', 'Nd', 'Nd', 'Nd', 'Nd', 'Nd', 'Nd', 'Nd', 'Ne', 'Ne', 'Ne', 'Ni', 'Ni', 'Ni', 'Ni', 'Ni', 'Np ', 'O', 'O', 'O', 'Os', 'Os', 'Os', 'Os', 'Os', 'Os', 'Os', 'P', 'Pb', 'Pb', 'Pb', 'Pb', 'Pd', 'Pd', 'Pd', 'Pd', 'Pd', 'Pd', 'Pr', 'Pt', 'Pt', 'Pt', 'Pt', 'Pt', 'Pt', 'Rb ', 'Rb', 'Re', 'Re', 'Rh', 'Ru', 'Ru', 'Ru', 'Ru', 'Ru', 'Ru', 'Ru', 'S', 'S', 'S', 'S', 'Sb', 'Sb', 'Sc', 'Se', 'Se', 'Se', 'Se', 'Se', 'Se', 'Si', 'Si', 'Si', 'Sm', 'Sm', 'Sm', 'Sm', 'Sm', 'Sm', 'Sm', 'Sn', 'Sn', 'Sn', 'Sn', 'Sn', 'Sn', 'Sn', 'Sn', 'Sn', 'Sn', 'Sr', 'Sr', 'Sr', 'Sr', 'Та', 'Та', 'Tb', 'Te', 'Te', 'Te', 'Te', 'Te', 'Te', 'Te', 'Te', 'Th', 'Ti', 'Ti', 'Ti', 'Ti', 'Ti', 'TI', 'TI', 'Tm', 'V', 'V', 'W', 'W', 'W', 'W', 'W', 'Yb', 'Yb', 'Yb', 'Yb', 'Yb', 'Yb', 'Yb', 'Zn', 'Zn', 'Zn', 'Zn', 'Zn', 'Zr', 'Zr', 'Zr', 'Zr', 'Zr ']

    #check for sulphur and nitrogen

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




    list_output = []
    mass_copy = mass.copy()
    abundance_copy = abundance.copy()
    isotopes_copy = isotopes.copy()

    for i in range (isotopes.count(list_atoms[0])):
        
        index = isotopes.index(list_atoms[0])
        list_output.append([mass_copy[index],abundance_copy[index]])
        mass_copy.pop(index)
        abundance_copy.pop(index)
        isotopes_copy.pop(index)

    list_atoms = list_atoms[1:]

    while len(list_atoms)>0:


        #This runs over all atoms in molecule

        list_output_new = []

        for i in range (len(list_output)):


            #This makes us run over all lists in list obtained before

            mass_copy = mass.copy()
            abundance_copy = abundance.copy()
            isotopes_copy = isotopes.copy() 

            for _ in range (isotopes.count(list_atoms[0])):


                #This for-loop runs over all isotope types of the atom type in pos 0 in list_atoms (input list)

                index = isotopes.index(list_atoms[0])
                new_mass = list_output[i][0] + mass_copy[index]
                new_proba = list_output[i][1] * abundance_copy[index]


                #removes any molecule who's probability is below 0.0001

                if render_imprecise_list: #only removes low-probability arrangements if render_imprecise_list arg is True
                    if new_proba>0.0001:
                        list_output_new.append([new_mass,new_proba])

                else:
                    list_output_new.append([new_mass,new_proba])

                mass_copy.pop(index)
                abundance_copy.pop(index)
                isotopes_copy.pop(index)
                
        list_output = list_output_new
        list_atoms.pop(0)



    #Conversion of list_output (which is a list of lists) to a combination of two lists (x_axis & y_axis)

    x_axis, y_axis = [],[]
    x_axis_final, y_axis_final = [],[]
    for j in range (len(list_output)):
        x_axis.append(round(list_output[j][0],3)) # Adds rounded value (should help with Python-limitations that render a diff of magnitude 10^(-7) to combinatorics
        
        y_axis.append(list_output[j][1]) #Adds the true value
        

    #Compression of lists x_axis & y_axis into x_axis_final & y_axis_final so that peaks corresponding to same mass will be represented together 
    
    for j in range (len(x_axis)):
            if x_axis.count(x_axis[j]) == 1:
                x_axis_final.append(x_axis[j])
                y_axis_final.append(y_axis[j])
            elif x_axis_final.count(x_axis[j]) == 0:
                x_axis_final.append(x_axis[j])
                y_axis_final.append(y_axis[j])
            else:
                index = x_axis_final.index(x_axis[j])
                y_axis_final[index] =y_axis_final[index] + y_axis[j]


    #if there is any, add peaks corresponding to Sulphur/Nitrogen presence
    maximum = max(y_axis_final)
    maximum_2 = 0
    for i in range (len(y_axis_final)):
        if y_axis_final[i]>maximum_2 and y_axis_final[i]< maximum:
            maximum_2 = y_axis_final[i]
    index = y_axis_final.index(maximum_2)

    if has_N:
        x_axis_final.append(x_axis_final[index] - 0.006)  
        y_axis_final.append(0.0035*count_N*maximum)  
    
    if has_S:
        x_axis_final.append(x_axis_final[index]-0.004)  
        y_axis_final.append(0.008*count_S*maximum)  

    return x_axis_final, y_axis_final



import unittest

class TestMainFunction(unittest.TestCase):

    def test_single_atom(self):
        list_atoms = ['H']
        imprecision_True_False = False
        result1, result2 = main_function(list_atoms, imprecision_True_False)
        expected_masses = [2.014, 1.008]
        expected_probabilities = [0.00015, 0.99985]
        self.assertEqual((result1, result2), (expected_masses, expected_probabilities))

    def test_single_atom_multiple(self):
        list_atoms = ['C', 'C']
        imprecision_True_False = False
        result1, result2 = main_function(list_atoms, imprecision_True_False)
        expected_masses = [26.007, 25.003, 24.0]
        expected_probabilities = [0.00012100000000000003, 0.021758000000000003, 0.9781210000000002]
        self.assertEqual((result1, result2), (expected_masses, expected_probabilities))

    def test_multiple_atoms(self):
        list_atoms = ['C', 'H', 'O']
        imprecision_True_False = False
        result1, result2 = main_function(list_atoms, imprecision_True_False)
        expected_masses = [33.017, 32.016, 31.012, 32.01, 31.01, 30.006, 32.013, 31.013, 30.009, 31.007, 30.007, 29.003]
        expected_probabilities =  [3.3e-09, 6.600000000000001e-10, 1.6460400000000003e-06, 2.1996700000000002e-05, 4.399340000000001e-06, 0.010971953960000001, 2.967e-07, 5.934e-08, 0.00014799396, 0.0019777033, 0.00039554066000000007, 0.9864784060400001]
        self.assertEqual((result1, result2), (expected_masses, expected_probabilities))

    def test_odd_count_nitrogen(self):
        list_atoms = ['N', 'N', 'N']
        imprecision_True_False = False
        result1, result2 = main_function(list_atoms, imprecision_True_False)
        expected_masses = [45.0, 44.003, 43.006, 42.009, 43.0]
        expected_probabilities = [5.0653000000000004e-08, 4.0918041e-05, 0.011018011958999999, 0.9889410193469999, 0.0103838807031435]
        self.assertEqual((result1, result2), (expected_masses, expected_probabilities))

    def test_even_count_nitrogen(self):
        list_atoms = ['N', 'N']
        imprecision_True_False = False
        result1, result2 = main_function(list_atoms, imprecision_True_False)
        expected_masses = [30.0, 29.003, 28.006]
        expected_probabilities = [1.3690000000000001e-05, 0.00737262, 0.9926136899999999]
        self.assertEqual((result1, result2), (expected_masses, expected_probabilities))

    def test_odd_count_sulfur(self):
        list_atoms = ['S', 'S', 'S']
        imprecision_True_False = False
        result1, result2 = main_function(list_atoms, imprecision_True_False)
        expected_masses = [107.901, 105.902, 104.906, 103.906, 103.903, 102.906, 101.907, 101.91, 100.911, 99.911, 101.904, 100.907, 99.908, 98.911, 97.912, 98.914, 97.915, 96.916, 95.916, 97.908]
        expected_probabilities = [8e-12, 5.052e-09, 9e-10, 1.1402399999999999e-07, 1.063446e-06, 3.789e-07, 4.800410400000001e-05, 3.375e-08, 8.5518e-06, 0.0005488323989999999, 7.4618461e-05, 3.9879225e-05, 0.005052431946, 0.0018001539, 0.11403374905199998, 4.21875e-07, 0.00016034625, 0.020314800899999996, 0.8579166140079998, 0.020589998736191994]
        self.assertEqual((result1, result2), (expected_masses, expected_probabilities))

    def test_imprecision_filter(self):
        list_atoms = ['C', 'C', 'H']
        imprecision_True_False = True
        result1, result2 = main_function(list_atoms, imprecision_True_False)
        expected_masses = [27.015, 26.011, 26.014, 25.008]
        expected_probabilities = [0.00012098185000000003, 0.021754736300000004, 0.00014671815000000002, 0.9779742818500002]
        self.assertEqual((result1, result2), (expected_masses, expected_probabilities))

if __name__ == '__main__':
    unittest.main()
