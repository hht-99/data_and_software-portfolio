#check for just prime factors
num_1=600851475143 #or int(input("Enter a number:"))
num_2=1

while num_2<=num_1:
    a = (num_1/num_2)
    num_2+=1
    if a - int(a) == 0:
        if a > 1:
            for j in range(2, int(a / 2) + 1):
                if (a % j) == 0:
                    break
            else:
                print(a, "is a prime factor")

#check for prime number
a = int(input("Enter a number:"))

if a > 1:
    for j in range(2, int(a/2) + 1):
        if (a % j) == 0:
            print(a, "is not a prime number")
            break
        else:
            print(a, "is a prime number")
else:
    print(a, "is not a prime number")

#check for prime and non-prime factors
num_1=int(input("Enter a number:"))
num_2=1

while num_2<=num_1:
    a = (num_1/num_2)
    num_2+=1
    if a - int(a) == 0:
        if a > 1:
            for j in range(2, int(a / 2) + 1):
                if (a % j) == 0:
                    print(a, "is a factor")
                    break
            else:
                print(a, "is a prime factor")
        else:
            print(a, "is a factor")