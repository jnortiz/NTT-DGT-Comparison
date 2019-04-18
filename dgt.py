#!/usr/env/python3
#coding: utf-8
import unittest
import time
from math import cos, sin, pi
from gaussian import GaussianInteger
from multiprocessing import Pool
from math import log

N = 2048

if N == 512:
    import n_512 as roots
else:
    if N == 1024:
        import n_1024 as roots
    else: # N = 2048
        import n_2048 as roots

p = 0xFFFFFFFF00000001 # 2**64 - 2**32 + 1
invMod = lambda y,q : pow(y,q-2,q)

nthroots = {
    256: GaussianInteger(4169533218321950981, 11340503865284752770),
    512: GaussianInteger(1237460157098848423, 590072184190675415),
    1024: GaussianInteger(13631489165933064639, 9250462654091849156),
}

invNthroots = {
    256: GaussianInteger(16440350318199496391, 8259263625103887671),
    512: GaussianInteger(11254465366323603399, 282547220712277683),
    1024: GaussianInteger(4772545667722300316, 8077569763565898552),
}

PROOTS = {
    0xFFFFFFFF00000001:7
}

######################################################################

def egcd(a, b):
    if a == 0:
        if type(b) == int:
            return (b, 0, 1)
        else:
            assert isinstance(b, GaussianInteger)
            return (b, GaussianInteger(0), GaussianInteger(1))
    else:
        g, x, y = egcd(b % a, a)
        return (g, y - (b // a) * x, x)

def modinv(b, n):
    g, x, _ = egcd(b, n)

    if g != 1:
        raise Exception('modular inverse does not exist (found %s)' % g)
    else:
        return x

def is_power2(n):
    n = int(n)
    while n>1:
        if n//2 != n//2.0: #compare integer division to float division
           return False
        n = n//2
    return True

######################
# Gentleman-sande DGT
#####################

# Apply DGT
def dgt_gentlemansande_online(x):    
    k = len(x) # n//2
    #assert is_power2(k)
    x = [a if isinstance(a, GaussianInteger) else GaussianInteger(a) for a in x]
    
    r = PROOTS[p] ## Primitive root of p

    #assert (p-1)%k == 0
    n = (p-1)//k

    g = pow(r, n, p)
    #assert pow(g, k, p) == 1 # k-th primitive root of unity
    gj = [pow(g, j, p) for j in range(k)]

    X = list(x) # Copy because the algorithm is run in-place
    for stride in range(int(log(k,2))):
        m = k // (2<<stride)

        for l in range(k // 2):
            j = l//(k//(2*m))
            a = pow(gj[j], k >> (int(log(k,2)) - stride), p)

            i = j + (l % (k//(2*m)))*2*m
        
            xi = X[i]
            xim = X[i + m]
            X[i] = xi + xim
            X[i + m] = a * (xi - xim)
    return X

def dgt_gentlemansande(x):    
    k = len(x) # n//2
    x = [a if isinstance(a, GaussianInteger) else GaussianInteger(a) for a in x]
    
    X = list(x) # Copy because the algorithm is run in-place
    for stride in range(int(log(k,2))):
        m = k // (2<<stride)

        for l in range(k // 2):
            j = l//(k//(2*m))
            a = pow(roots.gj[j], k >> (int(log(k,2)) - stride), p)

            i = j + (l % (k//(2*m)))*2*m
        
            xi = X[i]
            xim = X[i + m]
            X[i] = xi + xim
            X[i + m] = a * (xi - xim)
    return X    

def idgt_gentlemansande_online(x):
    k = len(x)
    #assert is_power2(k)
    x = [a if isinstance(a, GaussianInteger) else GaussianInteger(a) for a in x]

    r = PROOTS[p] ## Primitive root of p

    #assert (p-1)%k == 0
    n = (p-1)//k

    g = pow(r, n, p)
    #assert g != 0
    # print "g: %d" % g
    #assert pow(g, k, p) == 1 # n-th primitive root of unity
    invgj = [pow(g, (k - j), p) for j in range(k)] # g^-i \equiv g^((k-i) mod k) mod p

    X = list(x) # Copy because the algorithm is run in-place
    m = 1
    for stride in range(int(log(k,2))):
        for l in range(k // 2):
            j = l//(k//(2*m))
            #a = pow(n_512.gj_512_inv[j], k >> (stride + 1), p)
            a = pow(invgj[j], k >> (stride + 1), p)
            i = j + (l % (k//(2*m)))*2*m

            xi = X[i]
            xim = X[i + m]

            X[i] = xi + a * xim
            X[i + m] = xi - a * xim
        m = 2 * m
    return [v*modinv(k,p) for v in X]

def idgt_gentlemansande(x):
    k = len(x)
    x = [a if isinstance(a, GaussianInteger) else GaussianInteger(a) for a in x]

    X = list(x) # Copy because the algorithm is run in-place
    m = 1
    for stride in range(int(log(k,2))):
        for l in range(k // 2):
            j = l//(k//(2*m))
            a = pow(roots.invgj[j], k >> (stride + 1), p)
            i = j + (l % (k//(2*m)))*2*m

            xi = X[i]
            xim = X[i + m]

            X[i] = xi + a * xim
            X[i + m] = xi - a * xim
        m = 2 * m
    return [v*modinv(k,p) for v in X]

######################################################################

# Multipĺication in DGT's domain
def dgt_gentlemansande_mul_online(a, b):
    #assert len(a) == len(b)
    N = len(a)

    # Initialize
    a_folded = [GaussianInteger(x, y) for x, y in zip(a[:N//2], a[N//2:])]
    b_folded = [GaussianInteger(x, y) for x, y in zip(b[:N//2], b[N//2:])]

    # Compute h
    #assert pow(nthroots[N//2], N//2) == GaussianInteger(0, 1)
    #assert nthroots[N//2] * invNthroots[N//2] == GaussianInteger(1)

    # Twist the folded signals
    a_h = [a_folded[j] * pow(nthroots[N//2], j) for j in range(N // 2)]
    b_h = [b_folded[j] * pow(nthroots[N//2], j) for j in range(N // 2)]

    # Compute n//2 DGT
    a_dgt = dgt_gentlemansande_online(a_h)
    #assert idgt_gentlemansande_online(a_dgt) == a_h
    b_dgt = dgt_gentlemansande_online(b_h)
    #assert idgt_gentlemansande_online(b_dgt) == b_h

    # Point-wise multiplication
    c_dgt = [x * y for x, y in zip(a_dgt, b_dgt)]

    # Compute n//2 IDGT
    c_h = idgt_gentlemansande_online(c_dgt)

    # Remove twisting factors
    c_folded = [c_h[j] * pow(invNthroots[N//2], j) for j in range(N // 2)]

    # Unfold output
    c = [c_folded[j].re for j in range(N//2)] + [c_folded[j].imag for j in range(N//2)]

    return c

def dgt_gentlemansande_mul(a, b):
    
    N = len(a)

    # Initialize
    a_folded = [GaussianInteger(x, y) for x, y in zip(a[:N//2], a[N//2:])]
    b_folded = [GaussianInteger(x, y) for x, y in zip(b[:N//2], b[N//2:])]

    # Compute h
    #assert pow(nthroots[N//2], N//2) == GaussianInteger(0, 1)
    #assert nthroots[N//2] * invNthroots[N//2] == GaussianInteger(1)

    # Twist the folded signals
    a_h = [a_folded[j] * pow(nthroots[N//2], j) for j in range(N // 2)]
    b_h = [b_folded[j] * pow(nthroots[N//2], j) for j in range(N // 2)]

    # Compute n//2 DGT
    a_dgt = dgt_gentlemansande(a_h)
    #assert idgt_gentlemansande(a_dgt) == a_h
    b_dgt = dgt_gentlemansande(b_h)
    #assert idgt_gentlemansande(b_dgt) == b_h

    # Point-wise multiplication
    c_dgt = [x * y for x, y in zip(a_dgt, b_dgt)]

    # Compute n//2 IDGT
    c_h = idgt_gentlemansande(c_dgt)

    # Remove twisting factors
    c_folded = [c_h[j] * pow(invNthroots[N//2], j) for j in range(N // 2)]

    # Unfold output
    c = [c_folded[j].re for j in range(N//2)] + [c_folded[j].imag for j in range(N//2)]

    return c


# Apply schoolbook polynomial multiplication and reduces by Z_q // <x^N + 1>
def mul(a, b):
    assert len(a) == len(b)
    N = len(a)
    c = [0]*N

    # Mul and reduce
    for i in range(N):
        for j in range(N):
            v = a[i]*b[j]*(-1)**(int((i+j)//float(N)))

            c[(i+j) % N] = (c[(i+j) % N] + v) % p

    return c

# Multpĺication inside DGT's domain by an integer
def dgt_gentlemansande_mulscalar_online(a, b):
    #assert type(b) == int
    N = len(a)

    # Initialize
    a_folded = [GaussianInteger(x, y) for x, y in zip(a[:N//2], a[N//2:])]

    # Compute h
    #assert pow(nthroots[N//2], N//2) == GaussianInteger(0, 1)
    #assert nthroots[N//2] * invNthroots[N//2] == GaussianInteger(1)

    # Twist the folded signals
    a_h = [a_folded[j] * pow(nthroots[N//2], j) for j in range(N // 2)]

    # Compute n//2 DGT
    a_dgt = dgt_gentlemansande_online(a_h)
    #assert idgt_gentlemansande_online(a_dgt) == a_h

    # Point-wise multiplication
    c_dgt = [x * (b % p) for x in a_dgt]

    # Compute n//2 IDGT
    c_h = idgt_gentlemansande_online(c_dgt)

    # Remove twisting factors
    c_folded = [c_h[j] * pow(invNthroots[N//2], j) for j in range(N // 2)]

    # Unfold output
    c = [c_folded[j].re for j in range(N//2)] + [c_folded[j].imag for j in range(N//2)]

    return c

def dgt_gentlemansande_mulscalar(a, b):
    #assert type(b) == int
    N = len(a)

    # Initialize
    a_folded = [GaussianInteger(x, y) for x, y in zip(a[:N//2], a[N//2:])]

    # Compute h
    #assert pow(nthroots[N//2], N//2) == GaussianInteger(0, 1)
    #assert nthroots[N//2] * invNthroots[N//2] == GaussianInteger(1)

    # Twist the folded signals
    a_h = [a_folded[j] * pow(nthroots[N//2], j) for j in range(N // 2)]

    # Compute n//2 DGT
    a_dgt = dgt_gentlemansande(a_h)
    #assert idgt_gentlemansande(a_dgt) == a_h

    # Point-wise multiplication
    c_dgt = [x * (b % p) for x in a_dgt]

    # Compute n//2 IDGT
    c_h = idgt_gentlemansande(c_dgt)

    # Remove twisting factors
    c_folded = [c_h[j] * pow(invNthroots[N//2], j) for j in range(N // 2)]

    # Unfold output
    c = [c_folded[j].re for j in range(N//2)] + [c_folded[j].imag for j in range(N//2)]

    return c    

# Apply schoolbook polynomial multiplication and reduces by Z_q // <x^N + 1>
def mulint(a, b):
    assert type(b) == int
    N = len(a)
    c = [x * b % p for x in a]
    return c

############################################################## 

class TestDGTGentlemansande(unittest.TestCase):

    def test_transformation_online(self):
        # Verifies if iDGT(DGT(x)) == x        
        print("Testing DGT Gentleman-Sande Online")
        x = [x for x in range(N)]
        # print "\n".join([str(y) for y in dgt_gentlemansande(x)])
        start_time = time.time()
        idgt_gentlemansande_online(dgt_gentlemansande_online(x))
        end_time = time.time()
        print("----------", end_time - start_time, "s. ----------")

        self.assertEqual(idgt_gentlemansande_online(dgt_gentlemansande_online(x)), x)

    def test_transformation(self):
        # Verifies if iDGT(DGT(x)) == x        
        print("Testing DGT Gentleman-Sande")
        x = [x for x in range(N)]
        # print "\n".join([str(y) for y in dgt_gentlemansande(x)])
        start_time = time.time()
        idgt_gentlemansande(dgt_gentlemansande(x))
        end_time = time.time()
        print("----------", end_time - start_time, "s. ----------")

        self.assertEqual(idgt_gentlemansande(dgt_gentlemansande(x)), x)

    def test_mul_online(self):
        # Verifies multiplication in DGT's domain
        print("Polynomial multiplication using DGT Gentleman-Sande Online")
        a = [x for x in range(N)]
        b = [x for x in range(N)]

        start_time = time.time()
        dgt_gentlemansande_mul_online(a, b)
        end_time = time.time()

        print("----------", end_time - start_time, "s. ----------")

        self.assertEqual(
            dgt_gentlemansande_mul_online(a, b),
            mul(a, b)
            )

    def test_mul(self):
        # Verifies multiplication in DGT's domain
        print("Polynomial multiplication using DGT Gentleman-Sande")
        a = [x for x in range(N)]
        b = [x for x in range(N)]

        start_time = time.time()
        dgt_gentlemansande_mul(a, b)
        end_time = time.time()

        print("----------", end_time - start_time, "s. ----------")

        self.assertEqual(
            dgt_gentlemansande_mul(a, b),
            mul(a, b)
            )        

    def test_mulint_online(self):
        # Verifies multiplication in DGT's domain
        print("Multiplication by scalar using DGT Gentleman-Sande Online")
        a = [x for x in range(N)]
        b = p//3

        start_time = time.time()
        dgt_gentlemansande_mulscalar_online(a, b),
        end_time = time.time()

        print("----------", end_time - start_time, "s. ----------")

        self.assertEqual(
            dgt_gentlemansande_mulscalar_online(a, b),
            mulint(a, b)
            )


    def test_mulint(self):
        # Verifies multiplication in DGT's domain
        print("Multiplication by scalar using DGT Gentleman-Sande")
        a = [x for x in range(N)]
        b = p//3

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