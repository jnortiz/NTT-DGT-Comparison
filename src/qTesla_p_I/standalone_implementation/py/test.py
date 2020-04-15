import unittest
from gaussian import GaussianInteger
from dgt import idgt_gentleman_sande, dgt_cooley_tukey, idgt_cooley_tukey, dgt_gentleman_sande
import random
from params import N, p, kth_root_of_i, inv_kth_root_of_i, invkmodp

NRUNS = 1

fold   = lambda a: [GaussianInteger(x, y) for x, y in zip(a[:N//2], a[N//2:])]
unfold = lambda a: [a[j].re for j in range(N//2)] + [a[j].imag for j in range(N//2)]

def gen_polynomial_modp(length):
    x = []

    for _ in range(length):
        #x.append(random.randrange(0,p))
        x.append(1)
    
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

def poly_mul_ct_then_gs(x, y):

    # Folding the input signal into N//2 Gaussian integers
    x_folded = fold(x)
    # Applying the right-angle convolution
    x_folded = [x_folded[i] * kth_root_of_i[i] % p for i in range(N//2)]

    # Folding the input signal into N//2 Gaussian integers
    y_folded = fold(y)
    # Applying the right-angle convolution
    y_folded = [y_folded[i] * kth_root_of_i[i] % p for i in range(N//2)]

    x_dgt = dgt_cooley_tukey(x_folded)
    y_dgt = dgt_cooley_tukey(y_folded)

    xy_dgt = [(xi * yi) % p for xi, yi in zip(x_dgt, y_dgt)]

    xy = idgt_gentleman_sande(xy_dgt)

    # Removing the twisting factors and scaling by k^-1
    xy = [(xi * invkmodp * inv_kth_root_of_i[i]) % p for i, xi in enumerate(xy)]

    # Unfolding the output signal
    xy = unfold(xy)

    return xy

def poly_mul_gs_then_ct(x, y):

    # Folding the input signal into N//2 Gaussian integers
    x_folded = fold(x)
    # Applying the right-angle convolution
    x_folded = [x_folded[i] * kth_root_of_i[i] % p for i in range(N//2)]

    # Folding the input signal into N//2 Gaussian integers
    y_folded = fold(y)
    # Applying the right-angle convolution
    y_folded = [y_folded[i] * kth_root_of_i[i] % p for i in range(N//2)]

    x_dgt = dgt_gentleman_sande(x_folded)
    y_dgt = dgt_gentleman_sande(y_folded)

    xy_dgt = [xi * yi % p for xi, yi in zip(x_dgt, y_dgt)]

    xy = idgt_cooley_tukey(xy_dgt)

    # Removing the twisting factors and scaling by k^-1
    xy = [(xi * invkmodp * inv_kth_root_of_i[i]) % p for i, xi in enumerate(xy)]

    # Unfolding the output signal
    xy = unfold(xy)

    return xy

class TestDGT(unittest.TestCase):

    def _test_dgt_gs_then_ct(self):
        print("\nDGT transformation (Gentleman-Sande then Cooley-Tukey)")

        for _ in range(NRUNS):
            
            x = gen_polynomial_modp(N)

            # Folding the input signal into N//2 Gaussian integers
            x_folded = fold(x)
            # Applying the right-angle convolution
            x_folded_twisted = [x_folded[i] * kth_root_of_i[i] % p for i in range(N//2)]

            aux = dgt_gentleman_sande(x_folded_twisted)
            y_folded_twisted = idgt_cooley_tukey(aux)

            # Removing the twisting factors and scaling by k^-1
            y_folded = [(xi * invkmodp * inv_kth_root_of_i[i]) % p for i, xi in enumerate(y_folded_twisted)]

            # Unfolding the output signal
            y = unfold(y_folded)
           
            self.assertEqual(x, y)

    def _test_dgt_ct_then_gs(self):
        print("\nDGT transformation (Cooley-Tukey then Gentleman-Sande)")

        for _ in range(NRUNS):
            
            x = gen_polynomial_modp(N)

            # Folding the input signal into N//2 Gaussian integers
            x_folded = fold(x)
            # Applying the right-angle convolution
            x_folded_twisted = [x_folded[i] * kth_root_of_i[i] % p for i in range(N//2)]

            y_folded_twisted = idgt_gentleman_sande(dgt_cooley_tukey(x_folded_twisted))

            # Removing the twisting factors and scaling by k^-1
            y_folded = [(xi * invkmodp * inv_kth_root_of_i[i]) % p for i, xi in enumerate(y_folded_twisted)]

            # Unfolding the output signal
            y = unfold(y_folded)
           
            self.assertEqual(x, y)

    def _test_mul_gs_then_ct(self):
        print("\nPolynomial multiplication using DGT (Gentleman-Sande then Cooley-Tukey)")

        for _ in range(NRUNS):
            
            a = gen_polynomial_modp(N)
            b = gen_polynomial_modp(N)

            ab = schoolbook_mul(a, b)
            c = poly_mul_gs_then_ct(a, b)

            self.assertEqual(c, ab)

    def test_mul_ct_then_gs(self):
        print("\nPolynomial multiplication using DGT (Cooley-Tukey then Gentleman-Sande)")

        for _ in range(NRUNS):
            
            a = gen_polynomial_modp(N)
            b = gen_polynomial_modp(N)

            ab = schoolbook_mul(a, b)

            c = poly_mul_ct_then_gs(a, b)

            self.assertEqual(c, ab)

    def _test_basic_behaviour(self):

        x = [1]*N

        # Folding the input signal into N//2 Gaussian integers
        x_folded = fold(x)
        # Applying the right-angle convolution
        x_folded_twisted = [x_folded[i] * kth_root_of_i[i] % p for i in range(N//2)]

        aux = dgt_cooley_tukey(x_folded_twisted)

        print(aux)

        y_folded_twisted = idgt_gentleman_sande(aux)

        # Removing the twisting factors and scaling by k^-1
        y_folded = [(xi * invkmodp * inv_kth_root_of_i[i]) % p for i, xi in enumerate(y_folded_twisted)]

        # Unfolding the output signal
        y = unfold(y_folded)
        
        self.assertEqual(x, y)        


if __name__ == '__main__':
    unittest.main()