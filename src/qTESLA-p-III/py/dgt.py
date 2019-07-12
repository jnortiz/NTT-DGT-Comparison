import unittest
import time
from random import randint
from math import log

from gauss import GaussianInteger
from params import *
from params_dgt import *

## In-place DGT via Gentleman-Sande
def dgt_gentlemansande(x):   
    k = len(x)  
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

    return [v*invofkmodp for v in X]

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
    c_dgt = [(x * y) for x,y in zip(a_dgt, b_dgt)]

    # Compute n//2 IDGT
    c_h = idgt_gentlemansande(c_dgt)

    # Remove twisting factors
    c_folded = [c_h[j] * invNthroots[j] for j in range(N // 2)]

    # Unfold output
    c = [c_folded[j].re for j in range(N//2)] + [c_folded[j].imag for j in range(N//2)]

    return c

## Schoolbook polynomial multiplication in Z_q[x]/<x^N + 1>
def schoolbook_mul(a, b):
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
    a_h = [a_folded[j] * nthroots[j] for j in range(N // 2)]

    # Compute n//2 DGT
    a_dgt = dgt_gentlemansande(a_h)

    # Point-wise multiplication
    c_dgt = [x * (b % p) for x in a_dgt]

    # Compute n//2 IDGT
    c_h = idgt_gentlemansande(c_dgt)

    # Remove twisting factors
    c_folded = [c_h[j] * invNthroots[j] for j in range(N // 2)]

    # Unfold output
    c = [c_folded[j].re for j in range(N//2)] + [c_folded[j].imag for j in range(N//2)]

    return c

## Schoolbook scalar multiplication in Z_q[x]/<x^N + 1>
def mulint(a, b):
    assert type(b) == int
    N = len(a)
    c = [x * b % p for x in a]
    return c

def gen_polynomial_modp(length):
    x = []
    for i in range(length):
        #x.append(randint(0,p))
        x.append(1)
    return x

class TestDGTGentlemansande(unittest.TestCase):

    def test_transformation(self):
        print("\nTesting DGT Gentleman-Sande")        
        x = gen_polynomial_modp(N//2)        
        x = [a if isinstance(a, GaussianInteger) else GaussianInteger(a) for a in x]
        
        start_time = time.time()
        y = idgt_gentlemansande(dgt_gentlemansande(x))
        end_time = time.time()

        print(y)

        print("----------", end_time - start_time, "s. ----------")
        
        self.assertEqual(x, y)      

    def test_mul(self):
        print("\nPolynomial multiplication using DGT Gentleman-Sande")
        
        a = gen_polynomial_modp(N)
        b = gen_polynomial_modp(N)

        start_time = time.time()
        c = dgt_gentlemansande_mul(a, b)
        end_time = time.time()

        print("----------", end_time - start_time, "s. ----------")

        self.assertEqual(c, schoolbook_mul(a,b))

    def test_mulint(self):
        print("\nMultiplication by scalar using DGT Gentleman-Sande")
        
        a = gen_polynomial_modp(N)
        b = randint(0, p) # An arbitrary scalar number

        start_time = time.time()
        c = dgt_gentlemansande_mulscalar(a, b)
        end_time = time.time()

        print("----------", end_time - start_time, "s. ----------")

        self.assertEqual(c, mulint(a,b))

if __name__ == '__main__':
    unittest.main()