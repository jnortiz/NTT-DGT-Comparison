#!/usr/env/python3
#coding: utf-8
import sys
import unittest
import time
import random
from math import log,floor

parameter_set = 5

if parameter_set == 1:
    import qtesla_i as roots
    N = 512
    q = 4205569
else:
    if parameter_set == 2:
        import qtesla_iii_speed as roots
        N = 1024
        q = 8404993
    else:
        if parameter_set == 3:
            import qtesla_iii_size as roots
            N = 1024
            q = 4206593
        else:
            if parameter_set == 4:
                import qtesla_p_i as roots
                N = 1024
                q = 485978113
            else:
                if parameter_set == 5:
                    import qtesla_p_iii as roots
                    N = 2048
                    q = 1129725953
                else:
                    raise Exception("Please, choose a valid parameter set.")

invMod = lambda y,q:pow(y,q-2,q)

def ntt_cooley_tukey(x):
    a = list(x)
    t = N
    m = 1
    while m < N:
        t >>= 1
        for i in range(m):
            j1 = 2*i*t
            j2 = j1+t-1
            S = roots.psi_rev[m+i]
            for j in range(j1,j2+1):
                U = a[j]
                V = a[j+t]*S
                a[j] = (U+V)%q
                a[j+t] = (U-V)%q
        m<<=1
    return a

def intt_gentleman_sande(x):
    a = list(x)
    t = 1
    m = N
    while m > 1:
        j1 = 0
        h = m >> 1
        for i in range(h):
            j2 = j1 + t - 1
            S = roots.psi_inv_rev[h+i]
            for j in range(j1, j2+1):
                U = a[j]
                V = a[j+t]
                a[j] = (U+V)%q
                a[j+t] = ((U-V)*S)%q
            j1 = j1 + 2*t
        t <<= 1
        m >>= 1

    for j in range(N):
        a[j] = (a[j]*invMod(N,q))%q
    return a

def mul(a, b):
    assert len(a) == len(b)
    
    N = len(a)
    c = [0]*N

    for i in range(N):
        for j in range(N):
            v = a[i]*b[j]*(-1)**(int((i+j)//float(N)))
            c[(i+j) % N] = (c[(i+j) % N] + v) % q
    return c    

def NTT_Mul(a, b):
    a_ntt = ntt_cooley_tukey(a)
    b_ntt = ntt_cooley_tukey(b)

    for i in range(N):
        c_ntt = [(x*y)%q for x,y in zip(a_ntt, b_ntt)]

    c = intt_gentleman_sande(c_ntt) 
    return c   

def gen_polynomial_modq():
    x = []
    for i in range(N):
        x.append(random.randrange(0,q))
    return x

def test_polynomial_multiplication():
    a = gen_polynomial_modq()
    b = gen_polynomial_modq()
    c = NTT_Mul(a, b)
    return c

class TestNTT(unittest.TestCase):

    def test_transformation(self):
        a = []
        for i in range(N):
            a.append(random.randrange(0,q))
        assert(intt_gentleman_sande(ntt_cooley_tukey(a)) == a)

    def test_mul(self):        
        print("\nPolynomial multiplication using NTT Gentleman-Sande")

        a = []
        for i in range(N):
            a.append(random.randrange(0,q))

        b = []
        for i in range(N):
            b.append(random.randrange(0,q))

        start_time = time.time()        
        c = NTT_Mul(a, b)
        end_time = time.time()
        
        print("----------", end_time - start_time, "s. ----------")

        assert(mul(a,b) == c)

if __name__ == '__main__':
    unittest.main()