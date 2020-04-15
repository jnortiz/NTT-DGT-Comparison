import unittest
import time
from random import randint
from math import log

from gauss import GaussianInteger
from params import *
from params_dgt import *

modinv = lambda y,p:pow(y,p-2,p)

def dgt(x):
    n = len(x)
    x = [a if isinstance(a, GaussianInteger) else GaussianInteger(a) for a in x]

    g_inv = GaussianInteger(252694792,105937650)
    
    X = []
    for k in range(n):
        X.append(
            sum(
                [x[j]*pow(g_inv, j*k) % p for j in range(n)]
                )
            )
    return [x % p for x in X]

def idgt(X):
    n = len(X)
    invN = modinv(n, p)

    #assert (p-1)%n == 0
    #k = (p-1)//n

    #g = pow(r, k, p)
    #assert pow(g, n, p) == 1
    g = GaussianInteger(191706000,265394796)
    
    x = []
    for k in range(n):
        x.append(
            invN*sum(
                [X[j]*pow(g, j*k) % p for j in range(n)]
                ) % p
            )
    return x

def dgt_gentlemansande_mul(a, b):        
    N = len(a)

    # Initialize
    a_folded = [GaussianInteger(x, y) for x, y in zip(a[:N//2], a[N//2:])]
    b_folded = [GaussianInteger(x, y) for x, y in zip(b[:N//2], b[N//2:])]

    # Twist the folded signals
    a_h = [a_folded[j] * nthroots[j] for j in range(N // 2)] 
    b_h = [b_folded[j] * nthroots[j] for j in range(N // 2)]

    # Compute n//2 DGT
    a_dgt = dgt(a_h)
    b_dgt = dgt(b_h)

    # Point-wise multiplication
    c_dgt = [(x * y) for x,y in zip(a_dgt, b_dgt)]

    # Compute n//2 IDGT
    c_h = idgt(c_dgt)

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
    a_dgt = dgt(a_h)

    # Point-wise multiplication
    c_dgt = [x * (b % p) for x in a_dgt]

    # Compute n//2 IDGT
    c_h = idgt(c_dgt)

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
        x.append(randint(0,p))
    return x

class TestDGTGentlemansande(unittest.TestCase):

    def test_transformation(self):
        print("\nTesting DGT Gentleman-Sande")        
        x = gen_polynomial_modp(N//2)        
        x = [a if isinstance(a, GaussianInteger) else GaussianInteger(a) for a in x]
        
        start_time = time.time()
        y = idgt(dgt(x))
        end_time = time.time()

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