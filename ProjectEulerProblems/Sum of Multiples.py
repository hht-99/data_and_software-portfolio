#calc sum of all multiples of five and 3 under 1000

num = int(input("Enter a number: "))

if num > 0:
    sum = 0
    while (num > 0):
        if (num%3==0)or(num%5==0):
            sum += num
            num -= 1
        else:
            num -= 1
    print("The sum is", sum)
