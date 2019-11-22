from gaussian import GaussianInteger
from random import randint

modinvp = lambda y,p:pow(y,p-2,p)

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

def nthroot(n, p): 
    # First, find the factorization of p using https://www.alpertron.com.ar/GAUSSIAN.HTM
    # Then, try different combinations of signals and positions, real and imaginary parts, and invf0Modf1 and invf1Modf0 
    # (usually, the algorithm only find of inverse; the other is the conjugate of the first) in order to find a proper root of i.
    # Define the values of f0 and f1 below and obtain the generators.
    
    f0 = GaussianInteger(2841, 556) # Factorization of p
    f1 = GaussianInteger(2841, -556)

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

p = 8380417
r = 10
n = 256
k = n//2

def gen_roots_dgt():

    kth_root = None
    r = None            

    # 1. Compute the k-th root of i
    while r != 1:
        kth_root = nthroot(k, p) % p
        r, _, _ = egcd(kth_root, GaussianInteger(p))
    # 2. And its inverse mod p
    inv_kth_root = modinv(kth_root, GaussianInteger(p)) % p
    assert kth_root * inv_kth_root % p == 1

    pow_kth_root = []
    inv_pow_kth_root = []

    m = k
    w = kth_root
    w_inv = inv_kth_root
    while m > 1: # log(k) iterations            
        # 3. Compute the powers h, h^2, h^4, h^8, ..., h^k           
        pow_kth_root.append(w)
        inv_pow_kth_root.append(w_inv)
        
        w = (w * w) % p # h_{N}^{2} = h_{N/2}
        w_inv = (w_inv * w_inv) % p # h_{N}^{-2} = h_{N/2}^{-1}
        
        m >>= 1

    # 4. Invert the list of powers to h^k, ..., h^4, h^2, h
    inv_pow_kth_root.reverse()

    print("Powers of h: ", end="")
    print(pow_kth_root)

    print("\nPowers of h^{-1}: ", end="")
    print(inv_pow_kth_root)

gen_roots_dgt()

####################################

# nth_root = None
# r = None            

# while r != 1:
#     nth_root = nthroot(k, p) % p
#     r, _, _ = egcd(nth_root, GaussianInteger(p))

# inv_kth_root = modinv(nth_root, GaussianInteger(p)) % p
# assert nth_root * inv_kth_root % p == 1
# assert pow(nth_root, k) % p == GaussianInteger(0, 1)

# print(nth_root)
# print("--------------\n")
# print(inv_kth_root)