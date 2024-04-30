list = [0,1,2,3,4,5,6,7,8,9,10]
new_list = []

remove = True

for i in range (len(list)):
    if list[i]>4 and remove:
        new_list.append(list[i])

print(new_list)