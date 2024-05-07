x = [1,1.2,2,2.3,4,4.1,4.2]
y = [0,0,0,0,0,0,0]


def peak_merger(x_in, y_in):
    x_out, y_out = [],[]
    while len(x_in)>1:
        if x_in[0]>x_in[1]-0.2:
            y_in[1] = y_in[0] + y_in[1]
            x_in[1] = (x_in[0] + x_in[1])/2
            x_in.pop(0)
            y_in.pop(0)
        else:
            x_out.append(x_in[0])
            y_out.append(y_in[0])
            x_in.pop(0)
            y_in.pop(0)
    x_out.append(x_in[0])
    y_out.append(y_in[0])

    return x_out, y_out

print(peak_merger(x,y))
