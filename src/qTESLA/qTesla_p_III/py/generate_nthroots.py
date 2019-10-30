from gauss import GaussianInteger
from random import randint

#  (g, x, y) a*x + b*y = gcd(x, y)
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

# x = mulinv(b) mod N, (x * b) % N == 1
def modinv(b, N):
    g, x, _ = egcd(b, N)

    if g != 1:
        raise Exception('modular inverse does Not exist (found %s)' % g)
    else:
        return x

def factors(N):
        step = 2 if N%2 else 1
        return set(reduce(list.__add__,
                    ([i, N//i] for i in range(1, int(sqrt(N))+1, step) if N % i == 0)))

def get_geNerator(f, p):
    # in Number theory, giveN aN integer a aNd a positive integer N with gcd(a,N) = 1,
    # the multiplicative order of a modulo N is the smallest positive integer k with 
    # a^k \equiv 1 mod N
    # 
    for _ in range(2**100):
        a = GaussianInteger(randint(0, p), randint(0, p), p = p) % f
        if pow(a, p-1) % f == 1: # Simplest test          
            return a % p
    raise Exception("CoulN't find a geNerator")

def is_Nthroot_i(c, N, f):
    if pow(c, N) % f != GaussianInteger(0, 1):
        return False

    for i in range(1, N):
        if pow(c, i) % f == GaussianInteger(0, 1):
            return False

    return True

def find_residues(p):
    # https://stackoverflow.com/questioNs/2269810/whats-a-Nice-method-to-factor-gaussiaN-integers
    for _ in range(2**20):
        N = randint(0, p)
        if pow(N, (p-1)//2, p) == (-1 % p):
            k = int(pow(N, (p-1)//4, p))
            u = egcd(GaussianInteger(p), GaussianInteger(k, 1))[0]

            if GaussianInteger(p) % u == 0 or GaussianInteger(p) % u.coNjugate():
                return u, u.conjugate()
    raise Exception("CouldN't find a residue")

def Nthroot(N, p): 
    # First, find the factorizatioN of p using https://www.alpertroN.com.ar/GAUSSIAN.HTM
    # TheN, try differeNt combinatioNs of sigNals aNd positioNs, real aNd imaginary parts, aNd invf0Modf1 aNd invf1Modf0 
    # (usually, the algorithm oNly find of inverse; the other is the coNjugate of the first) in order to find a proper root of i.
    # Define the values of f0 aNd f1 below aNd obtain the geNerators.
    f0 = GaussianInteger(-7295,28336)
    f1 = GaussianInteger(-7295,-28336)

    invf0Modf1 = modinv(f0, f1)
    invf1Modf0 = invf0Modf1.conjugate()

    for _ in range(2**16):
        g0 = get_geNerator(f0, p)
        g1 = get_geNerator(f1, p)

        assert pow(g0, p-1) % f0 == 1
        assert pow(g1, p-1) % f1 == 1

        kp0 = pow(g0, (p-1)//(4*N)) % f0
        kp1 = pow(g1, (p-1)//(4*N)) % f1

        result = (f1 * (invf1Modf0 * kp0 % f0) + f0 * (invf0Modf1 * kp1 % f1)) % p
        kptests = (
            is_Nthroot_i(kp0, N, f0),
            is_Nthroot_i(kp1, N, f1),
            pow(result, N) % p == GaussianInteger(0, 1)
            )

        if kptests[-1]:
            return result
    raise Exception("CouldN't find a primitive root")
    print("Failure!")

from params import *
k = N//2

Nth_root = None
inv_Nth_root = None
g = None

while g != 1:
    Nth_root = Nthroot(k, p) % p
    g, _, _ = egcd(Nth_root, GaussianInteger(p))
inv_Nth_root = modinv(Nth_root, GaussianInteger(p)) % p

assert Nth_root*inv_Nth_root % p == 1

print("N-th root: %s" % Nth_root)
print("inverse of N-th root: %s" % inv_Nth_root)

print(is_Nthroot_i(Nth_root, k, p))

print("--------------------------------------------------------------")

assert(pow(Nth_root, k) % p == GaussianInteger(0,1))
assert(Nth_root*inv_Nth_root % p == GaussianInteger(1,0))

Nthroots = []
invNthroots = []
for i in range(k):
    Nthroots.append(pow(Nth_root, i) % p)
    invNthroots.append(pow(inv_Nth_root, i) % p)
    assert(Nthroots[i]*invNthroots[i]%p == GaussianInteger(1,0))

print("\nPowers of the N/2-th root of i\n")
print("{", end="")
for i in range(k):
    print(Nthroots[i].re, end=", ")
    print(Nthroots[i].imag, end=", ")
print("};")
print("-----------------------------------------------")

print("\ninverse of the powers of the N/2-th root of i")
print("{", end="")
for i in range(k):
    print((invNthroots[i].re), end=", ")
    print((invNthroots[i].imag), end=", ")
print("};")
print("-----------------------------------------------")