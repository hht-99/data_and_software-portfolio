#check for even number

def evennumbercheck(n_1=int(input("Enter a number:"))):
    if (n_1 / 2) - int(n_1 / 2) ==0:
        print('Even')
    else:
        print('Odd')


#call evennumbercheck(input)