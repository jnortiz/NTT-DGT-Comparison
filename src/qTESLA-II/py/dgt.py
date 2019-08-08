import unittest
import time
from random import randint
from math import log

from gauss import GaussianInteger
from params import *
from params_dgt import *

def poly_mul(x, y):
    # n = 2^7*6 = 768

    x_ntt = ntt_6(dgt_2(x))
    y_ntt = ntt_6(dgt_2(y))

    mul = [(x*y) for x,y in zip(x_ntt, y_ntt)]

    z = idgt_2(intt_6(mul))

    return z

def dgt_2(x):

    output = []

    # For each block of 2**7 coefficients
    for i in range(6):
        
        block = []

        # Extract the i-th block
        for j in range(2**7):
            block.append(x[(2**7)*i+j])

        # Fold the signal
        block = [GaussianInteger(x, y) for x, y in zip(block[:N//2], block[N//2:])]

        # Twist the folded signal
        block = [block[j]*nthroots[j] for j in range(N//2)]

        # Compute the DGT for each of them
        block = dgt_gentlemansande(block)

        # Save dgt(block) in the output without unfolding
        for j in range(2**6):
            output.append(block[j].re)
            output.append(block[j].imag)            

        # Save dgt(block) in the output <<with>> unfolding
        # for j in range(2**6):
        #     output.append(block[j].re)

        # for j in range(2**6):
        #     output.append(block[j].imag)

    print(output)
    return output   

def idgt_2(x):    
    output = []

    # For each block of 2**7 coefficients
    for i in range(6):
        block = []

        # Extract the block
        for j in range(2**7):
            block.append(x[(2**7)*i+j])

        xy_dgt = [GaussianInteger(block[j], block[j+1]) for j in range(2**6)]

        # Compute its counterpart in the ring defined by x^{2**7}+1
        xy = idgt_gentlemansande(xy_dgt)

        # Remove twisting factors
        xy = [xy[j]*invNthroots[j] for j in range(N//2)]

        # Unfold the output and save the idgt(block) into its correspondent place in the output
        output.append([xy[j].re for j in range(N//2)] + [xy[j].imag for j in range(N//2)])

    return output 

def ntt_6(X):

    x = list(X)

    for j in range(N):        
        x_tilde = []
        
        for i in range(6):
            x_tilde.append(x[N*i+j])

        x_folded = [x_tilde[i]*theta_powers[i] for i in range(6)] # theta transformation

        s0 = x_folded[0]+x_folded[3]
        s1 = x_folded[1]+x_folded[4]
        s2 = x_folded[2]+x_folded[5]
        s3 = (s1-s2)*zeta

        x[N*0+j] = (s0 + s1 + s2) % p
        x[N*2+j] = (s0 - s2 + s3) % p
        x[N*4+j] = (s0 - s3 - s1) % p

        s0 = x_folded[0] + x_folded[3]*theta_powers[3]
        s1 = x_folded[1]*theta_powers[1] + x_folded[4]*theta_powers[4]
        s2 = x_folded[2]*theta_powers[2] + x_folded[5]*theta_powers[5]
        s3 = (s1 - s2)*zeta

        x[N*1+j] = (s0 + s1 + s2) % p
        x[N*3+j] = (s0 - s2 + s3) % p 
        x[N*5+j] = (s0 - s3 - s1) % p

    return x

def intt_6(X):

    x = list(X)

    for j in range(N):        
        x_tilde = []
        
        for i in range(6):
            x_tilde.append(x[N*i+j])

        # mu_death
        t0 = (x_tilde[2] - x_tilde[4])*zeta
        t1 = (x_tilde[3] - x_tilde[5])*zeta

        s0 = x_tilde[0] + x_tilde[2] + x_tilde[4]
        s1 = x_tilde[1] + x_tilde[3] + x_tilde[5]

        x[N*0+j] = s0 - s1*zeta_squared
        x[N*3+j] = s0 - s1

        s0 = x_tilde[0] - x_tilde[2] - t0
        s1 = x_tilde[1] - x_tilde[3] - t1

        x[N*1+j] = s0 - s1*theta_powers[5]
        x[N*4+j] = s0 - s1*theta_powers[8]

        s0 = x_tilde[0] - x_tilde[4] + t0
        s1 = x_tilde[1] - x_tilde[5] + t1

        x[N*2+j] = s0 - s1*theta_powers[4]
        x[N*5+j] = s0 - s1*theta_powers[7]

        #psi_death

        x[N*0+j] = (x[N*0+j]*one_minus_zeta_by_nine) % p 
        x[N*1+j] = (x[N*1+j]*one_minus_zeta_by_nine*theta_inv) % p
        x[N*2+j] = (x[N*2+j]*one_minus_zeta_by_nine*theta_squared_inv) % p
        x[N*3+j] = (x[N*3+j]*zeta_minus_one_by_nine*theta_powers[3]) % p
        x[N*4+j] = (x[N*4+j]*zeta_minus_one_by_nine*theta_powers[2]) % p
        x[N*5+j] = (x[N*5+j]*zeta_minus_one_by_nine*theta_powers[1]) % p

    return x

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
#        x.append(randint(0,p))
        x.append(1)
    return x

class TestDGTGentlemansande(unittest.TestCase):

    def _test_transformation(self):
        print("\nTesting DGT Gentleman-Sande")        
        x = gen_polynomial_modp(N//2)        
        x = [a if isinstance(a, GaussianInteger) else GaussianInteger(a) for a in x]
        
        start_time = time.time()
        y = idgt_gentlemansande(dgt_gentlemansande(x))
        end_time = time.time()

        print("----------", end_time - start_time, "s. ----------")
        
        self.assertEqual(x, y)      

    def _test_mul(self):
        print("\nPolynomial multiplication using DGT Gentleman-Sande")
        
        a = gen_polynomial_modp(N)
        b = gen_polynomial_modp(N)

        start_time = time.time()
        c = dgt_gentlemansande_mul(a, b)
        end_time = time.time()

        print("----------", end_time - start_time, "s. ----------")

        self.assertEqual(c, schoolbook_mul(a,b))

    def _test_mulint(self):
        print("\nMultiplication by scalar using DGT Gentleman-Sande")
        
        a = gen_polynomial_modp(N)
        b = randint(0, p) # An arbitrary scalar number

        start_time = time.time()
        c = dgt_gentlemansande_mulscalar(a, b)
        end_time = time.time()

        print("----------", end_time - start_time, "s. ----------")

        self.assertEqual(c, mulint(a,b))

    def test_dgt_ntt(self):
        print("\nPolynomial multiplication in a non-power of two ring")
        
        a = gen_polynomial_modp(N*6)
        b = gen_polynomial_modp(N*6)

        start_time = time.time()
        c = poly_mul(a, b)
        end_time = time.time()

        print("----------", end_time - start_time, "s. ----------")

        #print(schoolbook_mul(a,b))

        print("\n\n")

        #print(c)

        self.assertEqual(c, schoolbook_mul(a,b))


if __name__ == '__main__':    
    unittest.main()