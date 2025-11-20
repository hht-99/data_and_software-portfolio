# largest palindrome made from the product of two 3-digit numbers

def backwards(x):
    return x[::-1]

for num_1 in range(100, 1000):
    for num_2 in range(100, 1000):
        product = num_1 * num_2
        check = str(product)
        test = backwards(check)
        issue = int(test)
        if product == issue and product>900000:
            print (product)