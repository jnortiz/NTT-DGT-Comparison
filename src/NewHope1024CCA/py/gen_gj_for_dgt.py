from math import log

#
# Generate g using Sage
# g = GF(p)
# g.zeta(k)
#
# Obtain g_inv sending a query to Wolfram Alpha
# inverse of g mod p
#

p = 12289
g = 3400
g_inv = 2859

N = 1024
k = N//2

def gen_powers_of_gj():    
    
    print("{", end="")
    for i in range(k//2):
        print((2**18*pow(g,i,p)) % p, end=", ")
    print("};")

def gen_powers_of_invgj():

    print("{", end="")
    for i in range(k//2):
        print((2**18*pow(g_inv,i,p)) % p, end=", ")
    print("};")

# def gen_powers_of_gj():    
#     gj = []
#     for i in range(k):
#         gj.append(pow(g,i,p))

#     print("[", end = ''),
#     for j in range(k-1): # For each gj[j]
#         print("[", end = ''),
#         for stride in range(int(log(k,2))-1): # Compute its powers mod p
#             a = pow(gj[j], k >> (int(log(k,2)) - stride), p)
#             print(a, end=", ")
#         a = pow(gj[j], k >> 1, p)
#         print(a, end="], ")
#         # print(a, end = 'u')
   
#     print("[", end =''),
#     for stride in range(int(log(k,2))-1): # Compute its powers mod p
#         a = pow(gj[k-1], k >> (int(log(k,2)) - stride), p)
#         print(a, end=", ")
#     a = pow(gj[k-1], k >> 1, p)
#     print(a, end="]]\n")

# def gen_powers_of_invgj():
#     invgj = []
#     for i in range(k):
#         invgj.append(pow(g_inv,i,p))

#     print("[", end = ''),
#     for j in range(k-1): # For each invgj[j]
#         print("[", end = ''),
#         for stride in range(int(log(k,2))-1): # Compute its powers mod p
#             a = pow(invgj[j], k >> (stride + 1), p)
#             print(a, end=", ")
#         a = pow(invgj[j], k >> (int(log(k,2))), p)
#         print(a, end="], ")
    
#     print("[", end = ''),
#     for stride in range(int(log(k,2))-1): # Compute its powers mod p
#         a = pow(invgj[k-1], k >> (stride + 1), p)
#         print(a, end=", ")
#     a = pow(invgj[k-1], k >> (int(log(k,2))), p)
#     print(a, end="]]\n")

gen_powers_of_gj()
print("\n\n-------------------------------------------------------\n\n")
gen_powers_of_invgj()