# Polynomial multiplication via DGT merged with twisting by powers of the k-th root of i mod p.

from gaussian import GaussianInteger
from params import N, p, invkmodp, w_inv, w
from random import randint
from math import log

fold   = lambda a: [GaussianInteger(x, y) for x, y in zip(a[:N//2], a[N//2:])]
unfold = lambda a: [a[j].re for j in range(N//2)] + [a[j].imag for j in range(N//2)]

def dgt_cooley_tukey(a):

    A = fold(a)
    n = len(A)

    m = 1
    Distance = n//2
    JTwiddle = 0
    
    while(m < n):
        
        for k in range(m):            
            JFirst = 2 * k * Distance
            JLast = JFirst + Distance - 1
            for j in range(JFirst, JLast+1, 1):
                Temp = w[JTwiddle] * A[j + Distance]
                A[j + Distance] = (A[j] - Temp) % p
                A[j] = (A[j] + Temp) % p
            JTwiddle += 1
        
        m = m * 2
        Distance = Distance//2

    return A

def idgt_gentleman_sande(a):

    A = list(a)
    n = len(A)

    ProblemSize = n
    Distance = 1

    h = 0
    while(ProblemSize > 1):
        
        for JFirst in range(Distance):
            j = JFirst
            JTwiddle = h
            while(j < n):
                Temp = A[j]
                A[j] = (Temp + A[j + Distance]) % p
                A[j + Distance] = (Temp - A[j + Distance]) * w_inv[JTwiddle] % p
                JTwiddle = JTwiddle + 1
                j = j + 2 * Distance

        h = h + ProblemSize//2
        ProblemSize = ProblemSize//2
        Distance = Distance * 2

    A = [ai * invkmodp % p for ai in A]

    return unfold(A)
