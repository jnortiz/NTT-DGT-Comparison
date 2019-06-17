 #!/usr/env/python3
#coding: utf-8
import unittest
import time
import random
from gauss import GaussianInteger
from math import log
from params import *
from params_dgt import *

## In-place DGT via Gentleman-Sande
def dgt_gentlemansande(x):   
    k = len(x) # n//2 | 256  
    X = list(x)
    
    for stride in range(int(log(k,2))):        
        m = k // (2<<stride)

        for l in range(k // 2):
            j = l//(k//(2*m))
            a = gj_powers[j][stride]
            i = j + (l % (k//(2*m)))*2*m
        
            xi = X[i]
            xim = X[i + m]
            X[i] = xi + xim
            X[i + m] = a * (xi - xim)

    return X    

## In-place IDGT via Cooley-Tukey
def idgt_gentlemansande(x):
    k = len(x)
    X = list(x)
    
    m = 1
    for stride in range(int(log(k,2))):
        for l in range(k // 2):
            j = l//(k//(2*m))            
            a = invgj_powers[j][stride]            
            i = j + (l % (k//(2*m)))*2*m
            xi = X[i]
            xim = X[i + m]

            X[i] = xi + a * xim
            X[i + m] = xi - a * xim
        m = 2 * m
    inv = 18374686475393433601

    return [v*inv for v in X]

def dgt_gentlemansande_mul(a, b):    
    
    N = len(a)

    # Initialize
    a_folded = [GaussianInteger(x, y) for x, y in zip(a[:N//2], a[N//2:])]
    b_folded = [GaussianInteger(x, y) for x, y in zip(b[:N//2], b[N//2:])]

    # Twist the folded signals
    a_h = [a_folded[j] * nthroots[j] for j in range(N // 2)] 
    b_h = [b_folded[j] * nthroots[j] for j in range(N // 2)]

    # Compute n//2 DGT
    a_dgt = dgt_gentlemansande(a_h)
    b_dgt = dgt_gentlemansande(b_h)

    # Point-wise multiplication
    c_dgt = [(x * y) for x, y in zip(a_dgt, b_dgt)]

    # Compute n//2 IDGT
    c_h = idgt_gentlemansande(c_dgt)

    # Remove twisting factors
    c_folded = [c_h[j] * invNthroots[j] for j in range(N // 2)]

    # Unfold output
    c = [c_folded[j].re for j in range(N//2)] + [c_folded[j].imag for j in range(N//2)]

    return c

## Schoolbook polynomial multiplication in Z_q[x]/<x^N + 1>
def mul(a, b):
    assert len(a) == len(b)
    N = len(a)
    c = [0]*N

    for i in range(N):
        for j in range(N):
            teste = (-1)**(int((i+j)//float(N)))
            v = a[i]*b[j]*teste
            c[(i+j) % N] = (c[(i+j) % N] + v) % p

    return c

## Scalar multipĺication in DGT's domain
def dgt_gentlemansande_mulscalar(a, b):
    N = len(a)

    # Initialize
    a_folded = [GaussianInteger(x, y) for x, y in zip(a[:N//2], a[N//2:])]

    # Twist the folded signals
    a_h = [a_folded[j] * nthroots[j] for j in range(N//2)]

    # Compute n//2 DGT
    a_dgt = dgt_gentlemansande(a_h)

    # Point-wise multiplication
    c_dgt = [x * (b % p) for x in a_dgt]

    # Compute n//2 IDGT
    c_h = idgt_gentlemansande(c_dgt)

    # Remove twisting factors
    c_folded = [c_h[j] * invNthroots[j] for j in range(N//2)]

    # Unfold output
    c = [c_folded[j].re for j in range(N//2)] + [c_folded[j].imag for j in range(N//2)]

    return c

## Schoolbook scalar multiplication in Z_q[x]/<x^N + 1>
def mulint(a, b):
    assert type(b) == int
    N = len(a)
    c = [x * b % p for x in a]
    return c

def gen_polynomial_modq(length):
    x = []
    for i in range(length):
        x.append(random.randrange(0,q))
    return x

def test_polynomial_multiplication():
    a = gen_polynomial_modq(N)
    b = gen_polynomial_modq(N)
    c = dgt_gentlemansande_mul(a, b)
    return c

def gen_powers_of_gj():
    k = len(gj)

    print("{", end = ''),
    for j in range(k-1): # For each gj[j]
        print("{", end = ''),
        for stride in range(int(log(k,2))-1): # Compute its powers mod p
            a = pow(gj[j], k >> (int(log(k,2)) - stride), p)
            print(a, end="u, ")
        a = pow(gj[j], k >> 1, p)
        print(a, end = 'u')
        print("},"),
    
    print("{", end = ''),
    for stride in range(int(log(k,2))-1): # Compute its powers mod p
        a = pow(gj[k-1], k >> (int(log(k,2)) - stride), p)
        print(a, end="u, ")
    a = pow(gj[k-1], k >> 1, p)
    print(a, end = 'u')
    print("}};")

def gen_powers_of_invgj():
    k = len(invgj)

    print("{", end = ''),
    for j in range(k-1): # For each invgj[j]
        print("{", end = ''),
        for stride in range(int(log(k,2))-1): # Compute its powers mod p
            a = pow(invgj[j], k >> (stride + 1), p)
            print(a, end="u, ")
        a = pow(invgj[j], k >> (int(log(k,2))), p)
        print(a, end="u")
        print("},"),
    
    print("{", end = ''),
    for stride in range(int(log(k,2))-1): # Compute its powers mod p
        a = pow(invgj[k-1], k >> (stride + 1), p)
        print(a, end="u, ")
    a = pow(invgj[k-1], k >> (int(log(k,2))), p)
    print(a, end="u")
    print("}};")    

class TestDGTGentlemansande(unittest.TestCase):

    def test_transformation(self):
        print("\nTesting DGT Gentleman-Sande")
        
        x = gen_polynomial_modq(N//2)        
        x = [a if isinstance(a, GaussianInteger) else GaussianInteger(a) for a in x]
        
        start_time = time.time()
        y = idgt_gentlemansande(dgt_gentlemansande(x))
        end_time = time.time()
        print("----------", end_time - start_time, "s. ----------")

        self.assertEqual(
            y,x)

    def test_mul(self):
        print("\nPolynomial multiplication using DGT Gentleman-Sande")
        
        a = gen_polynomial_modq(N)
        b = gen_polynomial_modq(N)

        start_time = time.time()
        c = dgt_gentlemansande_mul(a, b)
        end_time = time.time()

        print("----------", end_time - start_time, "s. ----------")

        self.assertEqual(
            dgt_gentlemansande_mul(a, b), mul(a,b)
            )

    def test_mulint(self):
        print("\nMultiplication by scalar using DGT Gentleman-Sande")
        
        a = gen_polynomial_modq(N)
        b = p//3 # An arbitrary scalar number

        start_time = time.time()
        dgt_gentlemansande_mulscalar(a, b),
        end_time = time.time()

        print("----------", end_time - start_time, "s. ----------")

        self.assertEqual(
            dgt_gentlemansande_mulscalar(a, b),
            mulint(a, b)
            )     

if __name__ == '__main__':
    unittest.main()
