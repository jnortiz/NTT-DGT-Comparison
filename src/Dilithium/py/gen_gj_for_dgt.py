from math import log

#
# Generate g using Sage
# g = GF(p)
# g.zeta(k)
#
# Obtain g_inv sending a query to Wolfram Alpha
# inverse of g mod p
#

p = 8380417
g = 3241972
g_inv = 6621070

N = 256
k = N//2

assert(pow(g,k,p) == 1)

def gen_powers_of_gj():    
    
    print("{", end="")
    for i in range(k):
        print(pow(g,i,p), end=", ")
    print("};")

def gen_powers_of_invgj():

    print("{", end="")
    for i in range(k):
        print(pow(g_inv,i,p), end=", ")
    print("};")

gen_powers_of_gj()
print("\n\n-------------------------------------------------------\n\n")
gen_powers_of_invgj()