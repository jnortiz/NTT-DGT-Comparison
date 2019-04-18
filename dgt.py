#!/usr/env/python3
#coding: utf-8
import unittest
from math import cos, sin, pi
from gaussian import GaussianInteger
from multiprocessing import Pool
from math import log
import time
import generate_prime as Prime
import n_512 as n_512

p = 0xFFFFFFFF00000001 # 2**64 - 2**32 + 1
invMod = lambda y,q : pow(y,q-2,q)
N = 512

nthroots = {
    4: GaussianInteger(17872535116924321793, 18446744035054843905),
    8: GaussianInteger(18446741870391328801, 18293621682117541889),
    16: GaussianInteger(13835058050988244993, 68719460336),
    32: GaussianInteger(2711615452560896026, 1702263320480097066),
    64: GaussianInteger(5006145799600815370, 13182758589097755684),
    128: GaussianInteger(1139268601170473317, 15214299942841048380),
    256: GaussianInteger(4169533218321950981, 11340503865284752770),
    512: GaussianInteger(1237460157098848423, 590072184190675415),
    1024: GaussianInteger(13631489165933064639, 9250462654091849156),
    2048: GaussianInteger(12452373509881738698, 10493048841432036102),
    4096: GaussianInteger(12694354791372265231, 372075061476485181),
    8192: GaussianInteger(9535633482534078963, 8920239275564743446),
    16384: GaussianInteger(9868966171728904500, 6566969707398269596),
    32768: GaussianInteger(10574165493955537594, 3637150097037633813),
    65536: GaussianInteger(2132094355196445245, 12930307277289315455)
}

invNthroots = {
    4: GaussianInteger(34359736320, 17868031517297999873),
    8: GaussianInteger(18311636080627023873, 18446741870391328737),
    16: GaussianInteger(18446739675663041537, 18446462594505048065),
    32: GaussianInteger(9223372049739972605, 9223372049739382781),
    64: GaussianInteger(3985917792403544159, 10871216858344511621),
    128: GaussianInteger(697250266387245748, 7269985899340929883),
    256: GaussianInteger(16440350318199496391, 8259263625103887671),
    512: GaussianInteger(11254465366323603399, 282547220712277683),
    1024: GaussianInteger(4772545667722300316, 8077569763565898552),
    2048: GaussianInteger(13028894352332048345, 9995848711590186809),
    4096: GaussianInteger(11525613835860693, 17335883825168514904),
    8192: GaussianInteger(17414606149056687587, 3916527805974289959),
    16384: GaussianInteger(9801605401783064476, 2749242888134484347),
    32768: GaussianInteger(10469048769509713349, 8715957816394874924),
    65536: GaussianInteger(15132804493885713016, 7997468840100395569)
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
def dgt_gentlemansande(x):    
    k = len(x) # n//2
    assert is_power2(k)
    x = [a if isinstance(a, GaussianInteger) else GaussianInteger(a) for a in x]
    
    r = PROOTS[p] ## Primitive root of p

    assert (p-1)%k == 0
    n = (p-1)//k

    g = pow(r, n, p)
    assert pow(g, k, p) == 1 # k-th primitive root of unity
    gj = [pow(g, j, p) for j in range(k)]

    X = list(x) # Copy because the algorithm is run in-place
    for stride in range(int(log(k,2))):
        m = k // (2<<stride)

        for l in range(k // 2):
            j = l//(k//(2*m))
            #a = pow(n_512.gj_512[j], k >> (int(log(k,2)) - stride), p)
            a = pow(gj[j], k >> (int(log(k,2)) - stride), p)

            i = j + (l % (k//(2*m)))*2*m
        
            xi = X[i]
            xim = X[i + m]
            X[i] = xi + xim
            X[i + m] = a * (xi - xim)
    return X

def idgt_gentlemansande(x):
    k = len(x)
    assert is_power2(k)
    x = [a if isinstance(a, GaussianInteger) else GaussianInteger(a) for a in x]

    r = PROOTS[p] ## Primitive root of p

    assert (p-1)%k == 0
    n = (p-1)//k

    g = pow(r, n, p)
    assert g != 0
    # print "g: %d" % g
    assert pow(g, k, p) == 1 # n-th primitive root of unity
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


######################################################################

# Multipĺication in DGT's domain
def dgt_gentlemansande_mul(a, b):
    assert len(a) == len(b)
    N = len(a)

    # Initialize
    a_folded = [GaussianInteger(x, y) for x, y in zip(a[:N//2], a[N//2:])]
    b_folded = [GaussianInteger(x, y) for x, y in zip(b[:N//2], b[N//2:])]

    # Compute h
    assert pow(nthroots[N//2], N//2) == GaussianInteger(0, 1)
    assert nthroots[N//2] * invNthroots[N//2] == GaussianInteger(1)

    # Twist the folded signals
    a_h = [a_folded[j] * pow(nthroots[N//2], j) for j in range(N // 2)]
    b_h = [b_folded[j] * pow(nthroots[N//2], j) for j in range(N // 2)]

    # Compute n//2 DGT
    a_dgt = dgt_gentlemansande(a_h)
    assert idgt_gentlemansande(a_dgt) == a_h
    b_dgt = dgt_gentlemansande(b_h)
    assert idgt_gentlemansande(b_dgt) == b_h

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
def dgt_gentlemansande_mulscalar(a, b):
    assert type(b) == int
    N = len(a)

    # Initialize
    a_folded = [GaussianInteger(x, y) for x, y in zip(a[:N//2], a[N//2:])]

    # Compute h
    assert pow(nthroots[N//2], N//2) == GaussianInteger(0, 1)
    assert nthroots[N//2] * invNthroots[N//2] == GaussianInteger(1)

    # Twist the folded signals
    a_h = [a_folded[j] * pow(nthroots[N//2], j) for j in range(N // 2)]

    # Compute n//2 DGT
    a_dgt = dgt_gentlemansande(a_h)
    assert idgt_gentlemansande(a_dgt) == a_h

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