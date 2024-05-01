
x = [1.2,1.4,1,1.8,2,2.4]
y = [2,2,1,2,2,2]




def list_sorter (x_in, y_in):
    x_out, y_out = [],[]
    print(len(x_in))
    while len(x_in)>0:
        min_x = min(x_in)
        index_min = x_in.index(min_x)
        x_out.append(min_x)
        y_out.append(y_in[index_min])
        x_in.pop(index_min)
        y_in.pop(index_min)
    print(len(x_out))
    return x_out, y_out

def peak_merger(x_in, y_in):
    x_out, y_out = [x_in[0]],[y_in[0]]

    for i in range(len(x_in)):    
        if x_in [i] > x_in[i-1]-0.3 and x_in[i]< x_in [i-1]+0.3:
            x_out.pop(i-1)
            x_out.append((x[i]+ x[i-1])/2)
            y_out.pop(i-1)
            y_out.append(y_in[i]+y_in[i-1])
        else:
            x_out.append(x_in[i])
            y_out.append(y_in[i])

    return x_out, y_out

print(peak_merger(x,y))




print(peak_merger(x,y))

