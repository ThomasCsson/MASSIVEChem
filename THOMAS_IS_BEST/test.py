list = [0.123,1.234,2.2585,3.2745,4.357875,5.35764,6,7,8,9,10]
new = []

for i in range(len(list)):
    new.append(round(list[i],2))

print(new)