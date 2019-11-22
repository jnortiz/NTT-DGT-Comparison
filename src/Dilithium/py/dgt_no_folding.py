# Work in progress. This version still does not work. It intendes to \
# merge the folding procedure with the polynomial multiplication.

from gaussian import GaussianInteger
from params import N, p, invkmodp, w_inv, w
from random import randint
from math import log

fold   = lambda a: [GaussianInteger(x, y) for x, y in zip(a[:N//2], a[N//2:])]
unfold = lambda a: [a[j].re for j in range(N//2)] + [a[j].imag for j in range(N//2)]

def dgt_cooley_tukey(a):

    A = [GaussianInteger(a[i], a[i+1]) for i in range(0, N, 2)]
    n = len(A)

    # NumOfGroups = n//2 = 64, for Dilithium
    JTwiddle = 63
    #Distance = 1
    for k in range(n//2):            
        JFirst = 2 * k
        JLast = JFirst
        for j in range(JFirst, JLast+1, 1):
            Temp = w[JTwiddle] * A[j+1]
            A[j+1] = (A[j] - Temp) % p
            A[j] = (A[j] + Temp) % p
        JTwiddle += 1

    NumOfGroups = 1
    Distance = n//2
    JTwiddle = 0    
    while(NumOfGroups < n and Distance >= 2):
        for k in range(NumOfGroups):            
            JFirst = 2 * k * Distance
            JLast = JFirst + Distance - 1
            for j in range(JFirst, JLast+1, 1):
                Temp = w[JTwiddle] * A[j + Distance]
                A[j + Distance] = (A[j] - Temp) % p
                A[j] = (A[j] + Temp) % p
            JTwiddle += 1
        
        NumOfGroups = NumOfGroups * 2
        Distance = Distance//2
    
    return A

def idgt_gentleman_sande(a):

    A = list(a)
    n = len(A)

    ProblemSize = n//2
    h = 64
    Distance = 2
    while(ProblemSize > 1):        
        for JFirst in range(Distance):
            j = JFirst
            JTwiddle = h        
            while(j < n):
                W = w_inv[JTwiddle]
                JTwiddle = JTwiddle + 1
                Temp = A[j]
                A[j] = (Temp + A[j + Distance]) % p
                A[j + Distance] = (Temp - A[j + Distance]) * W % p
                j = j + 2 * Distance

        h = h + ProblemSize//2
        ProblemSize = ProblemSize//2
        Distance = Distance * 2

    JTwiddle = 0
    j = 0
    while(j < n):
        W = w_inv[JTwiddle]
        JTwiddle += 1
        Temp = A[j]
        A[j] = (Temp + A[j + 1]) % p
        A[j + 1] = (Temp - A[j + 1]) * W % p
        j = j + 2
    
    A = [ai * invkmodp % p for ai in A]

    result = []
    for i in range(n):
        result.append(A[i].re)
        result.append(A[i].imag)
        
    return result

