x_axis = [10,11,12,12,13,13,13,14]
y_axis = [1,2,3,1,2,3,1,2]
x_axis_final, y_axis_final = [],[]



for j in range (len(x_axis)):
        if x_axis.count(x_axis[j]) == 1 or x_axis_final.count(x_axis[j]) == 0:
            x_axis_final.append(x_axis[j])
            y_axis_final.append(y_axis[j])
        else:
            index = x_axis_final.index(x_axis[j])
            y_axis_final[index] =y_axis_final[index] + y_axis[j]

print(x_axis_final,y_axis_final)