# Script for generating the powers of g multiplied by the powers of h for both forward and inverse DGT.
# The system parameters are those of Dilithium.

from gaussian import GaussianInteger

g = 3241972 # g^n = 1 mod p
g_inv = 6621070 # g * g_inv = 1 mod p
p = 8380417
N = 256
n = N//2 # Size of the transformation
R = 2**32 # Montgomery constant 

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

mth_root_of_i = {
    2:GaussianInteger(4616857, 3763560),
    4:GaussianInteger(3042578, 27908),
    8:GaussianInteger(4644566, 1062415),
    16:GaussianInteger(3120684, 3716310),
    32:GaussianInteger(2227391, 3639383),
    64:GaussianInteger(567084, 1577601),
    128:GaussianInteger(2719536, 8139413)
}

inv_mth_root_of_i = {
    128:GaussianInteger(3892641, 8377596),
    64:GaussianInteger(5638723, 2792435),
    32:GaussianInteger(3773247, 6178426),
    16:GaussianInteger(5491743, 7489321),
    8:GaussianInteger(6798825, 6106555),
    4:GaussianInteger(5206831, 158852),
    2:GaussianInteger(3763560, 3763560)
}

merged_wn_and_h = []
NumOfGroups = 1

while(NumOfGroups < n):    
    for k in range(NumOfGroups):
        merged_wn_and_h.append((brv_powers_of_g[k] * mth_root_of_i[2*NumOfGroups] * 2**32) % p)
    NumOfGroups = NumOfGroups * 2

#print(merged_wn_and_h)
#print("\n\n", len(merged_wn_and_h))

merged_wn_inv_and_h_inv = []

NumOfProblems = 1
ProblemSize = n

while(ProblemSize > 1):

    for JFirst in range(NumOfProblems):
        J = JFirst
        JTwiddle = 0

        while(J < n-1):
            W = (brv_powers_of_g_inv[JTwiddle] * inv_mth_root_of_i[ProblemSize] * R) % p
            merged_wn_inv_and_h_inv.append(W)
            print(W, end=", ")
            JTwiddle = JTwiddle + 1
            J = J + 2 * NumOfProblems
        print("\n")
    NumOfProblems = NumOfProblems * 2
    ProblemSize = ProblemSize//2
    print("\n")

#print(merged_wn_inv_and_h_inv)
#print(len(merged_wn_inv_and_h_inv))
