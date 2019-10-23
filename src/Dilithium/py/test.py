import unittest
from gaussian import GaussianInteger
from dgt import dgt, idgt
import random
from params import N, p

NRUNS = 1

def gen_polynomial_modp(length):
    x = []

    for i in range(length):
        #x.append(random.randrange(0,p))
        x.append(i)
    
    return x

def schoolbook_mul(a, b):
    assert len(a) == len(b)
    
    N = len(a)
    c = [0]*N

    for i in range(N):
        for j in range(N):
            v = a[i] * b[j] * (-1)**(int((i+j)//float(N)))
            c[(i + j) % N] = (c[(i + j) % N] + v) % p
    return c 

class TestDGTHierarquical(unittest.TestCase):

    def test_dgt(self):
        print("\nDGT transformation")

        for _ in range(NRUNS):
            x = gen_polynomial_modp(N)
            
            x_dgt = dgt(x)
            y = list(idgt(x_dgt))
            y = [y[i] % p for i in range(N)]
           
            self.assertEqual(x, y)

    def test_mul(self):
        print("\nPolynomial multiplication using DGT")

        for _ in range(NRUNS):
            #a = gen_polynomial_modp(N)
            #b = gen_polynomial_modp(N)
            a = [1]*N
            b = [1]*N

            ab = schoolbook_mul(a, b)

            print(ab)

            a_dgt = dgt(a)
            b_dgt = dgt(b)

            c_dgt = [x * y for x, y in zip(a_dgt, b_dgt)]

            c = idgt(c_dgt)

            self.assertEqual(c, ab)

if __name__ == '__main__':
    unittest.main()