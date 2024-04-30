
x = [1.2,1.4,1,1.8,2,2.4]
y = [2,2,2,2,2,2]

def peak_merger (x_axis_final, y_axis_final):
    x_out , y_out = [],[]
    for i in range (len(x_axis_final)):
        for j in range(i,len(x_axis_final)):
            if x_axis_final[i] < x_axis_final[j]+0.2 and x_axis_final[i] > x_axis_final[j]-0.2:
                x_out.append(x_axis_final[i])
                y_out.append(y_axis_final[i]+y_axis_final[j])
        
            
    return x_axis_final, y_axis_final


print(peak_merger(x,y))

