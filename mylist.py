mylist = [10, 'abc', 3, -5, 'Galiya']
print("Initial List:", mylist)

mylist.append(60)
print("After append(60):", mylist)



mylist.clear()
print("After clear():", mylist)


mylist.insert(2, 25)  
print("After insert(2, 25):", mylist)
mylist.insert(1,40)
mylist.insert(0,'abc')
print("after insert:",mylist)
mylist.remove(40)  
print("after removing 40:", mylist)
popped_element = mylist.pop(1)  
print("After pop(1):", mylist)
print("Popped Element:", popped_element)


my_list = [10, 20, 30, 40, 50]
print("\nReinitialized List:", my_list)



count_of_20 = my_list.count(20)
print("Count of 20:", count_of_20)


my_list.sort()  
print("After sort():", my_list)

my_list.reverse()
print("After reverse():", my_list)

copied_list = my_list.copy()  
print("Copied List:", copied_list)



length = len(my_list)  
print("Length of the list:", length)

