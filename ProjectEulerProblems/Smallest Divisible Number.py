# smallest positive number that is evenly divisible by all of the numbers from 1 to 20
# answer is 232792560

def simplecheck(num_a):
    n = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    for i in range(0, len(n)):
        if num_a % n[i] != 0:
            return 0
    return 1


for num_a in range(0, 1000000000):
    if simplecheck(num_a) == 1:
        print (num_a)