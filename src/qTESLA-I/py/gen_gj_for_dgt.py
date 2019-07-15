from math import log

p = 4205569
g = 2000313
g_inv = 730454

N = 512
k = N//2

def gen_powers_of_gj():        
    gj = []
    for j in range(k//2):
        gj.append(pow(g,j,p))

    for j in range(k//2):
        print(gj[j], end=", ")

def gen_powers_of_invgj():
    invgj = []
    for j in range(k//2):
        invgj.append(pow(g_inv,j,p))

    for j in reversed(range(k//2)):
        print(invgj[j], end=", ")

gen_powers_of_invgj()
