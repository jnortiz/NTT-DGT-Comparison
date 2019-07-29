from math import log

#
# Generate g using Sage
# g = GF(p)
# g.zeta(k)
#
# Obtain g_inv sending a query to Wolfram Alpha
# inverse of g mod p
#

p = 16801793
g = 8011249
g_inv = 4719840

N = 2048
k = N//2

def gen_powers_of_gj():    
    
    print("{", end="")
    for i in range(k//2):
        print(pow(g,i,p), end=", ")
    print("};")

def gen_powers_of_invgj():

    print("{", end="")
    for i in range(k//2):
        print(pow(g_inv,i,p), end=", ")
    print("};")

gen_powers_of_gj()
print("\n\n-------------------------------------------------------\n\n")
gen_powers_of_invgj()