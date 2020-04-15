from gaussian import GaussianInteger
from params import g, g_inv, p, N

n = N//2 # Size of the transformation

powers_of_g = []

for i in range(n//2):
    powers_of_g.append(pow(g, i, p))

powers_of_g_inv = []

for i in range(n//2):
    powers_of_g_inv.append(pow(g_inv, i, p))

def bit_reverse(x):
    X = list(x)
    N = len(X)
    j = 0
    for i in range(1, N):
        b = N >> 1
        while j >= b:
            j -= b
            b >>= 1
        j += b
        if j > i:
            X[i], X[j] = X[j], X[i]
    return X

brv_powers_of_g = bit_reverse(powers_of_g)
brv_powers_of_g_inv = bit_reverse(powers_of_g_inv)

brv_powers_of_g = [(2**32 * xi) % p for xi in brv_powers_of_g]
brv_powers_of_g_inv = [(2**32 * xi) % p for xi in brv_powers_of_g_inv]

# powers_of_g = [(2**32 * xi) % p for xi in powers_of_g]
# powers_of_g_inv = [(2**32 * xi) % p for xi in powers_of_g_inv]

print(powers_of_g)
print("\n-----------------------------------\n")
print(powers_of_g_inv)

# merged_wn_and_h = []
# NumOfGroups = 1

# while(NumOfGroups < n):    
#     for k in range(NumOfGroups):
#         merged_wn_and_h.append(brv_powers_of_g[k])
#     NumOfGroups = NumOfGroups * 2

# print(merged_wn_and_h)
# print("\n--------------------------------------\n")

# merged_wn_inv_and_h_inv = []

# NumOfProblems = 1
# ProblemSize = n

# while(ProblemSize > 1):

#     for JFirst in range(NumOfProblems):
#         J = JFirst
#         JTwiddle = 0

#         while(J < n-1):
#             W = brv_powers_of_g_inv[JTwiddle]
#             merged_wn_inv_and_h_inv.append(W)
#             print(W, end=", ")
#             JTwiddle = JTwiddle + 1
#             J = J + 2 * NumOfProblems
#         print("\n")
#     NumOfProblems = NumOfProblems * 2
#     ProblemSize = ProblemSize//2
#     print("\n")

# #print(merged_wn_inv_and_h_inv)