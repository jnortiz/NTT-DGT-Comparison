from math import log

#
# Generate g using Sage
# g = GF(p)
# g.zeta(k)
#
# Obtain g_inv sending a query to Wolfram Alpha
# inverse of g mod p
#

p = 8404993
g = 5384803
g_inv = 5920074

N = 1024
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
