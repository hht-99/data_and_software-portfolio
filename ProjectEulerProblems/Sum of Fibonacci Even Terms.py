#calc sum of all even terms in a fibonacci seq up to 4 mil

n_terms = int(32)

n_1 = 1
n_2 = 2
count = 0
sum = 0

print("The fibonacci sequence of the numbers is:")
while count < n_terms:
    print(n_1)
    nth = n_1 + n_2
    n_1 = n_2
    n_2 = nth
    count += 1
    if (n_1 / 2) - int(n_1 / 2) ==0:
        sum += n_1
print('The sum is', sum)