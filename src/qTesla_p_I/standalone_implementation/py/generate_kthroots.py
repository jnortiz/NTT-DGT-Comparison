from gaussian import GaussianInteger
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

# x = mulinv(b) mod n, (x * b) % n == 1
def modinv(b, n):
    g, x, _ = egcd(b, n)

    if g != 1:
        raise Exception('modular inverse does not exist (found %s)' % g)
    else:
        return x

def get_generator(f, p):
    # In number theory, given an integer a and a positive integer n with gcd(a,n) = 1,
    # the multiplicative order of a modulo n is the smallest positive integer k with 
    # a^k \equiv 1 mod n
    # 
    for _ in range(2**100):
        a = GaussianInteger(randint(0, p), randint(0, p), p = p) % f
        if pow(a, p-1) % f == 1: # Simplest test          
            return a % p
    raise Exception("Couln't find a generator")

def is_nthroot_i(c, n, f):
    if pow(c, n) % f != GaussianInteger(0, 1):
        return False

    for i in range(1, n):
        if pow(c, i) % f == GaussianInteger(0, 1):
            return False

    return True

def find_residues(p):
    # https://stackoverflow.com/questions/2269810/whats-a-nice-method-to-factor-gaussian-integers
    for _ in range(2**20):
        n = randint(0, p)
        if pow(n, (p-1)//2, p) == (-1 % p):
            k = int(pow(n, (p-1)//4, p))
            u = egcd(GaussianInteger(p), GaussianInteger(k, 1))[0]

            if GaussianInteger(p) % u == 0 or GaussianInteger(p) % u.conjugate():
                return u, u.conjugate()
    raise Exception("Couldn't find a residue")

def nthroot(n, p): 
    # First, find the factorization of p using https://www.alpertron.com.ar/GAUSSIAN.HTM
    # Then, try different combinations of signals and positions, real and imaginary parts, and invf0Modf1 and invf1Modf0 
    # (usually, the algorithm only find of inverse; the other is the conjugate of the first) in order to find a proper root of i.
    # Define the values of f0 and f1 below and obtain the generators.
    f0 = GaussianInteger(17281,6704) # Factorization of p
    f1 = GaussianInteger(17281,-6704)

    invf1Modf0 = modinv(f1, f0)
    invf0Modf1 = invf1Modf0.conjugate()

    for _ in range(2**16):
        g0 = get_generator(f0, p)
        g1 = get_generator(f1, p)

        assert pow(g0, p-1) % f0 == 1
        assert pow(g1, p-1) % f1 == 1

        kp0 = pow(g0, (p-1)//(4*n)) % f0
        kp1 = pow(g1, (p-1)//(4*n)) % f1

        result = (f1 * (invf1Modf0 * kp0 % f0) + f0 * (invf0Modf1 * kp1 % f1)) % p
        kptests = (
            is_nthroot_i(kp0, n, f0),
            is_nthroot_i(kp1, n, f1),
            pow(result, n) % p == GaussianInteger(0, 1)
            )

        if kptests[-1]:
            return result
    
    print("Failure!")
    raise Exception("Couldn't find a primitive root")

p = 343576577
n = 1024
k = n//2

nth_root = None
inv_nth_root = None
g = None

while g != 1:
    nth_root = nthroot(k, p) % p
    g, _, _ = egcd(nth_root, GaussianInteger(p))
inv_nth_root = modinv(nth_root, GaussianInteger(p)) % p

assert nth_root*inv_nth_root % p == 1

print("n-th root: %s" % nth_root)
print("inverse of n-th root: %s" % inv_nth_root)

print(is_nthroot_i(nth_root, k, p))

print("--------------------------------------------------------------")

assert(pow(nth_root, k) % p == GaussianInteger(0,1))
assert(nth_root*inv_nth_root % p == GaussianInteger(1,0))

nthroots = []
invNthroots = []
for i in range(k):
    nthroots.append(pow(nth_root, i) % p)
    invNthroots.append(pow(inv_nth_root, i) % p)
    assert(nthroots[i] * invNthroots[i] % p == GaussianInteger(1,0))

print("\nPowers of the n/2-th root of i\n")
print("{", end="")
for i in range(k):
    print((nthroots[i].re * 172048372) % p, end=", ")
    print((nthroots[i].imag * 172048372) % p, end=", ")
print("};")
print("-----------------------------------------------")

print("\nInverse of the powers of the n/2-th root of i")
print("{", end="")
for i in range(k):
    print((invNthroots[i].re * 172048372 * 342905529) % p, end=", ")
    print((invNthroots[i].imag * 172048372 * 342905529) % p, end=", ")
print("};")
print("-----------------------------------------------")