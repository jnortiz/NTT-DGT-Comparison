from gauss import GaussianInteger
import sys, traceback
from random import randint
from functools import reduce
from math import sqrt
import prime as Prime
import json
from fractions import gcd

#  (g, x, y) a*x + b*y = gcd(x, y)
def egcd(a, b):
    # print("a: %s(%d)\nb: %s(%d)\n" % (a, a.norm(), b, b.norm()))
    if a == 0:
        if type(b) == int:
            return (b, 0, 1)
        else:
            assert isinstance(b, GaussianInteger)
            return (b, GaussianInteger(0), GaussianInteger(1))
    else:
        g, x, y = egcd(b % a, a)
        # print(g.norm(), x.norm(), y.norm())
        return (g, y - (b // a) * x, x)

# x = mulinv(b) mod n, (x * b) % n == 1
def modinv(b, n):
    g, x, _ = egcd(b, n)

    if g != 1:
        raise Exception('modular inverse does not exist (found %s)' % g)
    else:
        return x

def factors(n):
        step = 2 if n%2 else 1
        return set(reduce(list.__add__,
                    ([i, n//i] for i in range(1, int(sqrt(n))+1, step) if n % i == 0)))

#############################################################$$$$$$$$$$$$$$$$$$
# nthroot
#

def get_generator(f, p):
    # In number theory, given an integer a and a positive integer n with gcd(a,n) = 1,
    # the multiplicative order of a modulo n is the smallest positive integer k with 
    # a^k \equiv 1 mod n
    # 
    for _ in range(2**100):
        a = GaussianInteger(randint(0, p), randint(0, p), p = p) % f
        if pow(a, p-1) % f == 1: # Simplest test
            # For doing it correctly we have to execute this second test
            # however, it is too computationally complex. 
            # Given the context we can trust the simplest test and return.
            # If the result is not trully a generator we handle it later!
            # 
            # if len(set(
            # [(x.re, x.imag) for x in [pow(a,i) for i in range(p)]]
            # )) == p - 1: 
            
            return a % p
    raise Exception("Couln't find a generator")

def is_nthroot_i(c, n, f):
    # Verifies if c is a nth-root of i mod f
    if pow(c, n) % f != GaussianInteger(0, 1):
        # print("is_nthroot_i: Failed 1")
        return False

    for i in range(1, n):
        if pow(c, i) % f == GaussianInteger(0, 1):
            # print("is_nthroot_i: Failed 2")
            return False

    return True

def find_residues(p):
    # 
    # https://stackoverflow.com/questions/2269810/whats-a-nice-method-to-factor-gaussian-integers
    #
    for _ in range(2**20):
        n = randint(0, p)
        if pow(n, (p-1)//2, p) == (-1 % p):
            k = int(pow(n, (p-1)//4, p))
            # if k % 2 == 1:
                # return GaussianInteger(k, 1), GaussianInteger(k, -1)
            u = egcd(GaussianInteger(p), GaussianInteger(k, 1))[0]

            if GaussianInteger(p) % u == 0 or GaussianInteger(p) % u.conjugate():
                return u, u.conjugate()
    raise Exception("Couldn't find a residue")

def nthroot(n, p): 
    # 
    # Let n bt a positive integer. A primitivie nth-rooth of i is an nth root
    # of i that is not a kth root of i for any positive k < n.
    # 
    # if p == (2**64 - 2**32 + 1):
    #     f0 = GaussianInteger(65536,  4294967295)
    #     f1 = GaussianInteger(65536, -4294967295)
    #     invf0Modf1 = GaussianInteger(0, 2147483648)
    #     invf1Modf0 = GaussianInteger(0,-2147483648)
    # elif p == 524269:
    #     f0 = GaussianInteger(262,  675, p)
    #     f1 = GaussianInteger(262, -675, p)
    #     invf0Modf1 = GaussianInteger(0, 321940, p)
    #     invf1Modf0 = GaussianInteger(0,-321940, p)
    # else:
        # assert p == 523777
        # f0 = GaussianInteger(284,  665, p)
        # f1 = GaussianInteger(284, -665, p)
        # invf0Modf1 = GaussianInteger(0, 224485, p)
        # invf1Modf0 = GaussianInteger(0,-224485, p)
    # assert p % 4 == 1
    # invf0Modf1 = None
    # invf1Modf0 = None
    # f0, f1 = find_residues(p)
    # try:
    #     invf0Modf1 = modinv(f0, f1)
    #     invf1Modf0 = modinv(f1, f0)
    #     assert f0 * f1 == p
    #     assert (f0 * invf0Modf1) % f1 == 1
    #     assert (f1 * invf1Modf0) % f0 == 1
    #     print("Found: f0, f1 = %s, %s" % (f0, f1))
    #     print("Found: invf0Modf1, invf1Modf0 = %s, %s" % (invf0Modf1, invf1Modf0))
    # except Exception as e:
    #     invf0Modf1 = None
    #     invf1Modf0 = None
    #     print("Bad residues.")
    # if f0 is None or f1 is None or invf0Modf1 is None or invf1Modf0 is None:
    #     raise Exception("Can't find invertible residues.")

    f0 = GaussianInteger(-1663, 1200)
    f1 = GaussianInteger(-1663, -1200)
    invf0Modf1 = GaussianInteger(-895, -512)
    invf1Modf0 = GaussianInteger(-895, 512)

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
        # print("kp0: %s, kp1: %s, kp: %s" % kptests)

        # If result == GaussianInteger(0, 1) then we found a nthroot of i
        # Otherwise, try another pair of generators 
        if kptests[-1]:
            return result
        # else:
            # print("Failure: %s != %s" % (pow(result, n), GaussianInteger(0, 1, p = p)))
    raise Exception("Couldn't find a primitive root")


    print("Failure!")

from fractions import gcd
import prime as Prime

# n = 63
arange, brange = 4, 15
quantity = 20
# Receives a integer and a list.
# Returns true if n is coprime with every integer in the list
def coprime_check(n, l):
    for i in l:
        if gcd(n,i) != 1:
            return False
    return True

def gi_check(x):
    if x % 2 == 1 and x % 4 == 1: # avoid 2 and integers non congruents to 1 mod 4
        for n in range(arange, brange):
            # avoid integers not divisible by n and 4*n
            if (x-1) % (4*pow(2, n)) != 0 or (x-1) % pow(2, n) != 0:
                return False
        return True

qi = 4205569
n = 256

nthroots = {}
invnthroots = {}

g = None
while g != 1:
    nthroots = nthroot(n, qi) % qi
    g, _, _ = egcd(nthroots, GaussianInteger(qi))
invnthroots = modinv(nthroots, GaussianInteger(qi)) % qi

assert nthroots*invnthroots % qi == 1

print("nthroots: %s" % nthroots)
print("invnthroots: %s" % invnthroots)

print("--------------------------------------------------------------")

p = 4205569
N = 512
k = N//2

nth_root = GaussianInteger(1299482,3370617)
inv_nth_root = GaussianInteger( 2720832,2389903)

assert(pow(nth_root, k) % p == GaussianInteger(0,1))
assert(nth_root*inv_nth_root % p == GaussianInteger(1,0))

nthroots = []
invNthroots = []
for i in range(k):
    nthroots.append(pow(nth_root, i) % p)
    invNthroots.append(pow(inv_nth_root, i) % p)
    assert(nthroots[i]*invNthroots[i]%p == GaussianInteger(1,0))
    print(nthroots[i])

print("\nPowers of the n/2-th root of i\n")
print("[", end="")
for i in range(k):
    print("GaussianInteger(", end="")
    print(nthroots[i].re, end=",")
    print(nthroots[i].imag, end="), ")
print("]")
print("-----------------------------------------------")

print("\nInverse of the powers of the n/2-th root of i")
print("[", end="")
for i in range(k):
    print("GaussianInteger(", end="")
    print(invNthroots[i].re, end=",")
    print(invNthroots[i].imag, end="), ")
print("]")
print("-----------------------------------------------")
