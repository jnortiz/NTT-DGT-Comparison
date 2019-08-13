from math import log

p = 33564673
g = 6766491
g_inv = 27777950

N = 18
k = N//2

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
    
def gen_powers_of_gj():        
    for j in range(k):
        print(pow(g,j,p), end=", ")

def gen_powers_of_invgj():
    for j in range(k):
        print(pow(g_inv,j,p), end=", ")


gen_powers_of_gj()
print("\n\n----------------------------------------\n\n")
gen_powers_of_invgj()
print("\n\n")
